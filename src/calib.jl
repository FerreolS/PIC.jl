using Parameters

@with_kw struct BboxParams
    BBOX_DX_LOWER::Int = 2
    BBOX_DX_UPPER::Int = 2
    BBOX_DY_LOWER::Int = 21
    BBOX_DY_UPPER::Int = 18
    BBOX_WIDTH::Int = BBOX_DX_LOWER + 1 + BBOX_DX_UPPER
    BBOX_HEIGHT::Int = BBOX_DY_LOWER + 1 + BBOX_DY_UPPER
end

@with_kw struct PICParams{R<:Real,Q}
    @deftype R
    nλ::Int = 3
    @assert (nλ == 3 || nλ == 4)
    NLENS::Int = 18908
    @assert NLENS ≥ 1
    lasers_order::Int = 2
    @assert lasers_order ≥ 1
    lasers_λs::Vector{R} = [987.72e-9, 1123.71e-9, 1309.37e-9, 1545.10e-9][1:nλ]
    LASERS_CXY0S_INIT_PATH::String = joinpath(dirname(pathof(PIC)), "lasers_cxy0s_init.txt")
    lasers_cxy0s_init::Matrix{R} = readdlm(LASERS_CXY0S_INIT_PATH, Float64)
    @assert size(lasers_cxy0s_init) == (NLENS, 2)
    @assert all(isfinite.(lasers_cxy0s_init))
    LASERS_CX1_INIT = -0.6001811340726275
    LASERS_CX2_INIT = -0.3187688427580339
    LASERS_CY1_INIT = 89.9795748752424
    LASERS_CY2_INIT = -52.635157560302524

    lasers_fwhms_init::Vector{R} = [2.3, 2.4, 2.7, 2.9][1:nλ]
    @assert length(lasers_fwhms_init) == nλ
    lamp_order::Int = 2

    λLAMP_RANGE::Q = LinRange(850e-9, 1600e-9, 10000) # coarse wavelength range of the instrument

    lamp_cfwhms_init::Vector{R} = [2.5, 0, 0, 0][1:(lamp_order+1)]

    bbox_params::BboxParams = BboxParams()
end


function fitSpectralLawAndProfile(
    lasers::WeightedArray,
    lamp::WeightedArray,
    ; calib_params::PICParams,
    valid_lenslets::AbstractVector{Bool}=trues(NLENS)
)
    @unpack_PICParams calib_params
    @unpack_BboxParams bbox_params

    size(valid_lenslets) == (NLENS,) || throw(ArgumentError("valid_lenslets must be of size NLENS"))

    λref = mean(lasers_λs)

    bboxes = fill(BoundingBox{Int}(-1, -1, -1, -1), NLENS)
    lasers_cxs = fill(NaN64, lasers_order + 1, NLENS)
    lasers_cys = fill(NaN64, lasers_order + 1, NLENS)
    lasers_fwhms = fill(NaN64, nλ, NLENS)
    lasers_amplitudes = fill(NaN64, nλ, NLENS)
    lasers_pixels_dists = fill(NaN64, BBOX_WIDTH, BBOX_HEIGHT, NLENS)
    lasers_pixels_λs = fill(NaN64, BBOX_WIDTH, BBOX_HEIGHT, NLENS)
    lamp_cfwhms = fill(NaN64, lamp_order + 1, NLENS)
    lamp_cxs = fill(NaN64, lamp_order + 1, NLENS)
    lamp_backs = fill(NaN64, NLENS)
    lamp_amplitudes = fill(NaN64, BBOX_HEIGHT, NLENS)

    p = Progress(NLENS; showspeed=true)

    assigned_lenslets = falses(NLENS)

    @inbounds for i in findall(valid_lenslets)

        bbox = get_bbox(lasers_cxy0s_init[i, 1], lasers_cxy0s_init[i, 2]; bbox_params=bbox_params)
        if !ismissing(bbox)
            bboxes[i] = bbox
            assigned_lenslets[i] = true
        end
    end

    Threads.@threads for i in findall(valid_lenslets)
        try

            # lasers

            lens_lasers = view(lasers, bboxes[i])

            lasers_cxs_init = [lasers_cxy0s_init[i, 1];
                LASERS_CX1_INIT * (λref * 1e6);
                LASERS_CX2_INIT * (λref * 1e6)^2]

            lasers_cys_init = [lasers_cxy0s_init[i, 2];
                LASERS_CY1_INIT * (λref * 1e6);
                LASERS_CY2_INIT * (λref * 1e6)^2]

            lasers_lkl = Lasers_LKL(
                nλ, lasers_order, lasers_λs, λref, bboxes[i], lens_lasers)

            (fit_lasers_cxs, fit_lasers_cys, fit_fwhms, fit_amplitudes) = fit_lens_lasers(lasers_lkl,
                lasers_fwhms_init, lasers_cxs_init, lasers_cys_init)

            lasers_cxs[:, i] .= fit_lasers_cxs
            lasers_cys[:, i] .= fit_lasers_cys
            lasers_fwhms[:, i] .= fit_fwhms
            lasers_amplitudes[:, i] .= fit_amplitudes

            lens_lasers_pixels_dists = view(lasers_pixels_dists, :, :, i)
            lens_lasers_pixels_λs = view(lasers_pixels_λs, :, :, i)



            compute_lasers_dists_and_λmap!(
                λLAMP_RANGE, bboxes[i], lasers_order, λref, fit_lasers_cxs, fit_lasers_cys,
                lens_lasers_pixels_dists, lens_lasers_pixels_λs)

            # lamp

            lens_lamp = view(lamp, bboxes[i])

            lamp_cxs_init = [fit_lasers_cxs[1]; 0; 0]

            lamp_lkl = Lamp_LKL(lamp_order, λref, bboxes[i], lens_lamp, lens_lasers_pixels_dists, lens_lasers_pixels_λs)

            (fit_lamp_cfwhms, fit_lamp_cxs, fit_lamp_back, fit_lamp_amplitudes) = fit_lens_lamp(
                lamp_lkl, lamp_cfwhms_init, lamp_cxs_init)

            lamp_cfwhms[:, i] .= fit_lamp_cfwhms
            lamp_cxs[:, i] .= fit_lamp_cxs
            lamp_backs[i] = fit_lamp_back
            lamp_amplitudes[:, i] .= fit_lamp_amplitudes

        catch e
            @debug "Error on lenslet $i" exception = (e, catch_backtrace())
            assigned_lenslets[i] = false
        end
        next!(p)
    end
    ProgressMeter.finish!(p)

    (; nλ, lasers_λs, λref, lasers_order, lamp_order, assigned_lenslets, bboxes,
        lasers_cxs, lasers_cys, lasers_fwhms, lasers_amplitudes,
        lasers_pixels_dists, lasers_pixels_λs, lamp_cfwhms, lamp_cxs, lamp_backs, lamp_amplitudes)
end

function get_bbox(center_x::Float64, center_y::Float64; bbox_params::BboxParams=BboxParams())
    @unpack_BboxParams bbox_params
    bbox = round(
        Int,
        BoundingBox(; xmin=center_x - BBOX_DX_LOWER,
            xmax=center_x + BBOX_DX_UPPER,
            ymin=center_y - BBOX_DY_LOWER,
            ymax=center_y + BBOX_DY_UPPER),
        RoundNearestTiesUp) # rounding mode to preserve bbox size

    size(bbox) == (BBOX_WIDTH, BBOX_HEIGHT) || return missing
    ((bbox.xmin ≥ 1) & (bbox.xmax ≤ 2048) & (bbox.ymin ≥ 1) & (bbox.ymax ≤ 2048)) || return missing
    return bbox
end

const NLENS = 18908

const LASERS_ORDER_DEFAULT = 2

const LASERS_4_λS = [987.72e-9, 1123.71e-9, 1309.37e-9, 1545.10e-9]
const LASERS_3_λS = LASERS_4_λS[1:3]

const LASERS_CXY0S_INIT_PATH = joinpath(dirname(pathof(PIC)), "lasers_cxy0s_init.txt")
const LASERS_CXY0S_INIT = readdlm(LASERS_CXY0S_INIT_PATH, Float64)
const LASERS_CX1_INIT = -0.6001811340726275
const LASERS_CX2_INIT = -0.3187688427580339
const LASERS_CY1_INIT = 89.9795748752424
const LASERS_CY2_INIT = -52.635157560302524

const LASERS_FWHMS_INIT = [2.3, 2.4, 2.7]

const LAMP_ORDER_DEFAULT = 2

const λLAMP_RANGE = LinRange(850e-9, 1600e-9, 10000) # coarse wavelength range of the instrument

const LAMP_CFWHMS_INIT = [2.3, 2.5, 2.9]

const BBOX_DX_LOWER = 2
const BBOX_DX_UPPER = 2
const BBOX_DY_LOWER = 21
const BBOX_DY_UPPER = 18
const BBOX_WIDTH = BBOX_DX_LOWER + 1 + BBOX_DX_UPPER
const BBOX_HEIGHT = BBOX_DY_LOWER + 1 + BBOX_DY_UPPER

function fitSpectralLawAndProfile(
    lasers::WeightedArray,
    lamp::WeightedArray,
    ; nλ::Int,
    lasers_fwhms_init::Vector{Float64},
    lasers_order::Int=LASERS_ORDER_DEFAULT,
    lasers_cxy0s_init::Matrix{Float64}=LASERS_CXY0S_INIT,
    lamp_order::Int=LAMP_ORDER_DEFAULT,
    lamp_cfwhms_init::Vector{Float64}=LAMP_CFWHMS_INIT,
    valid_lenslets::AbstractVector{Bool}=trues(NLENS)
)
    NLENS ≥ 1 || throw(ArgumentError)
    nλ ≥ 2 || throw(ArgumentError)
    size(lasers_fwhms_init) == (nλ,) || throw(ArgumentError)
    size(lasers_cxy0s_init) == (NLENS, 2) || throw(ArgumentError)
    lamp_order ≥ 1 || throw(ArgumentError)
    length(lamp_cfwhms_init) == lamp_order + 1 || throw(ArgumentError)
    size(valid_lenslets) == (NLENS,) || throw(ArgumentError)

    lasers_λs = (nλ == 4) ? LASERS_4_λS : (nλ == 3) ? LASERS_3_λS : throw(ArgumentError)
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

    assigned_lenslets = copy(valid_lenslets)

    Threads.@threads for i in findall(valid_lenslets)
        try

            bbox = get_bbox(lasers_cxy0s_init[i, 1], lasers_cxy0s_init[i, 2])
            bboxes[i] = bbox

            # lasers

            lens_lasers = view(lasers, bbox)

            lasers_cxs_init = [lasers_cxy0s_init[i, 1];
                LASERS_CX1_INIT * (λref * 1e6);
                LASERS_CX2_INIT * (λref * 1e6)^2]

            lasers_cys_init = [lasers_cxy0s_init[i, 2];
                LASERS_CY1_INIT * (λref * 1e6);
                LASERS_CY2_INIT * (λref * 1e6)^2]

            lasers_lkl = Lasers_LKL(
                nλ, lasers_order, lasers_λs, λref, bbox, lens_lasers)

            (fit_lasers_cxs, fit_lasers_cys, fit_fwhms, fit_amplitudes) = fit_lens_lasers(lasers_lkl,
                lasers_fwhms_init, lasers_cxs_init, lasers_cys_init)

            lasers_cxs[:, i] .= fit_lasers_cxs
            lasers_cys[:, i] .= fit_lasers_cys
            lasers_fwhms[:, i] .= fit_fwhms
            lasers_amplitudes[:, i] .= fit_amplitudes

            lens_lasers_pixels_dists = view(lasers_pixels_dists, :, :, i)
            lens_lasers_pixels_λs = view(lasers_pixels_λs, :, :, i)



            compute_lasers_dists_and_λmap!(
                λLAMP_RANGE, bbox, lasers_order, λref, fit_lasers_cxs, fit_lasers_cys,
                lens_lasers_pixels_dists, lens_lasers_pixels_λs)

            # lamp

            lens_lamp = view(lamp, bbox)

            lamp_cxs_init = [fit_lasers_cxs[1]; 0; 0]

            lamp_lkl = Lamp_LKL(lamp_order, λref, bbox, lens_lamp, lens_lasers_pixels_dists, lens_lasers_pixels_λs)

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

function get_bbox(center_x::Float64, center_y::Float64)::BoundingBox{Int}

    bbox = round(
        Int,
        BoundingBox(; xmin=center_x - BBOX_DX_LOWER,
            xmax=center_x + BBOX_DX_UPPER,
            ymin=center_y - BBOX_DY_LOWER,
            ymax=center_y + BBOX_DY_UPPER),
        RoundNearestTiesUp) # rounding mode to preserve bbox size

    size(bbox) == (BBOX_WIDTH, BBOX_HEIGHT) || error()
    ((bbox.xmin ≥ 1) & (bbox.xmax ≤ 2048) & (bbox.ymin ≥ 1) & (bbox.ymax ≤ 2048)) || error()

    bbox
end

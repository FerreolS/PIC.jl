using Parameters

@with_kw struct BboxParams
    BBOX_DX_LOWER::Int = 2
    BBOX_DX_UPPER::Int = 2
    BBOX_DY_LOWER::Int = 21
    BBOX_DY_UPPER::Int = 18
    BBOX_WIDTH::Int = BBOX_DX_LOWER + 1 + BBOX_DX_UPPER
    BBOX_HEIGHT::Int = BBOX_DY_LOWER + 1 + BBOX_DY_UPPER
end


@with_kw struct OptimParams{R<:Real,Q}
    @deftype R
    maxeval::Int = 500
    xtol::Tuple{R,R} = (0.0, 1e-7)
    ftol::Tuple{R,R} = (0.0, 1e-8)
    gtol::Tuple{R,R} = (0.0, 1e-6)
    verb::Bool = false
    lower::R = -Inf
    upper::R = +Inf
    ADbackend::Q = AutoZygote()
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

    laserOptim::OptimParams = OptimParams()
    lampOptim::OptimParams = OptimParams()

    multi_thread::Bool = true
end

struct LensletModel{T}
    λref::T
    bbox::BoundingBox{Int}
    lasers_cxs::Vector{T}
    lasers_cys::Vector{T}
    lasers_fwhms::Vector{T}
    λs::Vector{T}
    fwhm_coefs::Vector{T}
    x_coefs::Vector{T}
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

    bboxes = fill(BoundingBox{Int}(nothing), NLENS)
    lenslet_array = Vector{LensletModel{Float64}}(undef, NLENS)

    lasers_amplitudes = fill(NaN64, nλ, NLENS)
    lasers_pixels_dists = Vector{Vector{Float64}}(undef, NLENS)
    lamp_backs = fill(NaN64, NLENS)
    lamp_amplitudes = fill(NaN64, BBOX_HEIGHT, NLENS)
    lasers_cost = fill(NaN64, NLENS)
    lamp_cost = fill(NaN64, NLENS)
    lasers_model = zeros(size(lasers))
    lamp_model = zeros(size(lamp))
    p = Progress(NLENS; showspeed=true)

    assigned_lenslets = falses(NLENS)

    @inbounds for i in findall(valid_lenslets)

        bbox = get_bbox(lasers_cxy0s_init[i, 1], lasers_cxy0s_init[i, 2]; bbox_params=bbox_params)
        if !ismissing(bbox)
            bboxes[i] = bbox
            assigned_lenslets[i] = true
        end
    end

    #Threads.@threads for i in findall(assigned_lenslets)
    # from https://discourse.julialang.org/t/optionally-multi-threaded-for-loop/81902/8?u=skleinbo
    _foreach = multi_thread ? ThreadsX.foreach : Base.foreach
    _foreach(findall(assigned_lenslets)) do i
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

            (fit_lasers_cxs, fit_lasers_cys, fit_fwhms, fit_amplitudes, model, cost) = fit_lens_lasers(lasers_lkl,
                lasers_fwhms_init, lasers_cxs_init, lasers_cys_init; optim=laserOptim)

            lasers_cost[i] = cost
            view(lasers_model, bboxes[i]) .= model
            lasers_amplitudes[:, i] .= fit_amplitudes

            lens_lasers_pixels_dists = fill(NaN64, BBOX_HEIGHT)
            lens_lasers_pixels_λs = fill(NaN64, BBOX_HEIGHT)



            compute_lasers_λmap!(
                λLAMP_RANGE, bboxes[i], lasers_order, λref, fit_lasers_cxs, fit_lasers_cys,
                lens_lasers_pixels_dists, lens_lasers_pixels_λs)

            # lamp

            lens_lamp = view(lamp, bboxes[i])

            lamp_cxs_init = [fit_lasers_cxs[1]; 0; 0]

            lamp_lkl = Lamp_LKL(lamp_order, λref, bboxes[i], lens_lamp, lens_lasers_pixels_dists, lens_lasers_pixels_λs)

            (fit_lamp_cfwhms, fit_lamp_cxs, fit_lamp_back, fit_lamp_amplitudes, model, cost) = fit_lens_lamp(
                lamp_lkl, lamp_cfwhms_init, lamp_cxs_init)


            lenslet_array[i] = LensletModel{Float64}(λref,
                bboxes[i], fit_lasers_cxs, fit_lasers_cys, fit_fwhms,
                lens_lasers_pixels_λs, fit_lamp_cfwhms, fit_lamp_cxs)

            view(lamp_model, bboxes[i]) .= model
            lamp_cost[i] = cost


            lamp_backs[i] = fit_lamp_back
            lamp_amplitudes[:, i] .= fit_lamp_amplitudes

            lasers_pixels_dists[i] = lens_lasers_pixels_dists


        catch e
            @debug "Error on lenslet $i" exception = (e, catch_backtrace())
            assigned_lenslets[i] = false
        end
        next!(p)
    end
    ProgressMeter.finish!(p)

    (; lenslet_array, nλ, lasers_λs, λref, lasers_order, lamp_order, assigned_lenslets, bboxes,
        lasers_amplitudes,
        lasers_pixels_dists, lamp_backs, lamp_amplitudes,
        lasers_cost, lamp_cost, lasers_model, lamp_model)
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


function extract_spectrum(wd::WeightedArray{T,N},
    lenslet::LensletModel;
    restrict=0.01,
    nonnegative=false,
    kwds...
) where {T,N}
    bbox = lenslet.bbox
    if N > 2
        (; data, precision) = view(wd, bbox, :)
    else
        (; data, precision) = view(wd, bbox)
    end

    model = T.(compute_lamp_images(lenslet))

    model = model ./ sum(model, dims=1)
    if restrict > 0
        precision .*= (model .> restrict)
    end

    αprecision = sum(model .^ 2 .* precision, dims=1)
    α = sum(model .* precision .* data, dims=1) ./ αprecision

    nanpix = .!isnan.(α)
    if nonnegative
        positive = nanpix .& (α .>= T(0))
    else
        positive = nanpix
    end

    wp = WeightedArray(dropdims(positive .* α, dims=1), dropdims(positive .* αprecision, dims=1))

    return wp
end

function extract_spectra(wd::WeightedArray{T,N},
    lenslets::AbstractVector{<:LensletModel},
    assigned_lenslets::AbstractVector{Bool};
    multi_thread::Bool=true,
    restrict=0.01,
    nonnegative=false,
    kwds...) where {T,N}

    wspectra = Vector{WeightedArray{T}}(undef, length(lenslets))

    l = size(lenslets[findfirst(assigned_lenslets)].bbox, 2)
    unseen = WeightedArray(zeros(l), zeros(l))
    _foreach = multi_thread ? ThreadsX.foreach : Base.foreach
    _foreach(eachindex(lenslets, assigned_lenslets)) do i
        if assigned_lenslets[i]
            wspectra[i] = extract_spectrum(wd, lenslets[i]; restrict=restrict, nonnegative=nonnegative)
        else
            wspectra[i] = unseen
        end
    end
    return wspectra
end

function extract_wavelegth_coordinates(lenslet_array::AbstractVector{<:LensletModel}, assigned_lenslets::AbstractVector{Bool}; thrs=3)
    idx = findall(assigned_lenslets)
    isnothing(idx) && throw(ArgumentError("No assigned lenslets"))
    λs = fill(NaN64, size(lenslet_array[findfirst(assigned_lenslets)].bbox, 2), length(lenslet_array))
    map(i -> (λs[:, i] .= lenslet_array[i].λs), idx)
    # mapreduce(x -> x.λs, hcat, lenslet_array[assigned_lenslets])
    #assigned_lenslets .&= identify_bad_lenslet(λs)

    blue = λs[1, idx]
    red = λs[end, idx]
    Δλ = red .- blue
    view(assigned_lenslets, idx) .&= ((blue .- median(blue)) .< thrs * mad(blue)) .&&
                                     ((red .- median(red)) .< thrs * mad(red)) .&&
                                     ((Δλ .- median(Δλ)) .< thrs * mad(Δλ))

    return λs, assigned_lenslets
end


function build_λrange(λs::AbstractMatrix{<:Real}; superres=1)
    blue = λs[1, :]
    red = λs[end, :]
    return range(start=minimum(blue), stop=maximum(red), step=median((red .- blue) ./ 40) / superres)
end


function find_index(knots::AbstractRange, sample)
    return (sample - first(knots)) / step(knots) + 1
end

function build_sparse_interpolation_matrix(knots, samples; kernel::Kernel{T,N}=CatmullRomSpline()) where {T,N}
    lk = length(kernel)
    lin = length(samples)
    col = length(knots)

    nelement = lk * lin
    L = zeros(Int, nelement)
    C = zeros(Int, nelement)
    V = zeros(T, nelement)
    c = 1

    for (l, sample) ∈ enumerate(samples)
        offweights = InterpolationKernels.compute_offset_and_weights(kernel, T.(find_index(knots, sample)))
        weights = vcat(offweights[2]...)
        off::Int = round(Int, offweights[1]) + 1
        L[c:(c+lk-1)] .= l
        C[c:(c+lk-1)] .= min.(max.(off:(off+lk-1), 1), col)
        V[c:(c+lk-1)] .= weights
        c += lk
    end
    return sparse(L, C, V, lin, col)
end

function get_lower_uppersamples(λ::AbstractVector)
    lower = [(3 * λ[1] .- λ[2]) / 2; (λ[1:end-1] .+ λ[2:end]) / 2]
    upper = [(λ[2:end] .+ λ[1:end-1]) / 2; (3 * λ[end] .- λ[end-1]) / 2]
    return lower, upper
end


function reverse_cumsum(v)
    out = similar(v)
    out[end] = v[end]
    @inbounds for i ∈ (length(v)-1):-1:1
        out[i] = out[i+1] + v[i]
    end
    return out
end

function build_sparse_interpolation_integration_matrix(knots, lowersample, uppersamples; kernel::Kernel{T,N}=CatmullRomSpline()) where {T,N}

    lk = length(kernel)
    lin = length(uppersamples)
    lin == length(lowersample) || throw(DimensionMismatch("uppersamples and lowersample must have the same length"))
    col = length(knots)

    nelement = col * lin
    L = zeros(Int, nelement)
    C = zeros(Int, nelement)
    V = zeros(T, nelement)
    c = 1

    for (l, (lsample, usample)) ∈ enumerate(zip(lowersample, uppersamples))
        uoffweights = InterpolationKernels.compute_offset_and_weights(kernel, T.(find_index(knots, usample)))
        loffweights = InterpolationKernels.compute_offset_and_weights(kernel, T.(find_index(knots, lsample)))
        uweights = vcat(uoffweights[2]...)[2:end]
        uoff::Int = round(Int, uoffweights[1]) + 1

        lweights = vcat(loffweights[2]...)[2:end]
        loff::Int = round(Int, loffweights[1]) + 1

        lv = uoff - loff + lk - 1
        v = ones(T, lv)
        v[(lv-lk+2):end] .= reverse_cumsum(uweights)
        v[1:lk-1] .-= reverse_cumsum(lweights)
        off = min.(max.(loff+1:(loff+lv), 1), col)
        L[c:(c+lv-1)] .= l
        C[c:(c+lv-1)] .= off
        V[c:(c+lv-1)] .= v
        c += lv
    end
    return sparse(L[1:c-1], C[1:c-1], V[1:c-1], lin, col)
end
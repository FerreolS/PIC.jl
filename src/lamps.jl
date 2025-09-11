
@concrete struct Lamp_LKL
    order::Int
    λref <: Real   # reference wavelength
    bbox::BoundingBox{Int}
    data <: WeightedArray
    lens_lasers_pixels_dists <: AbstractVector{<:Real}
    lasers_pixels_λs <: AbstractVector{<:Real}
end




function fit_lens_lamp(
    lamp_lkl::Lamp_LKL,
    lamp_cfwhms_init::Vector{Float64},
    lamp_cxs_init::Vector{Float64};
    optim=OptimParams()
)
    @unpack_OptimParams optim

    vmlmbvars = encode_lamp_lkl_vmlmbvars(lamp_cfwhms_init, lamp_cxs_init)
    grad = similar(vmlmbvars)
    prep = prepare_gradient(lamp_lkl, ADbackend, vmlmbvars)
    fg!(x, grad) = DifferentiationInterface.value_and_gradient!(lamp_lkl, grad, prep, ADbackend, x)[1]
    vmlmb!(fg!, vmlmbvars; verb=verb, maxeval=maxeval, ftol=ftol, xtol=xtol, gtol=gtol, lower=lower, upper=upper)

    #vmlmb!(lamp_lkl, vmlmbvars; verb=false, ftol=(0.0, 1e-8), maxeval=500, autodiff=true)

    (fit_cfwhms, fit_cxs) = decode_lamp_lkl_vmlmbvars(vmlmbvars)

    (cost, fit_back, fit_amplitudes, model) = compute_lamp_cost_and_back_and_amplitudes(
        lamp_lkl, fit_cfwhms, fit_cxs)

    (copy(fit_cfwhms), copy(fit_cxs), fit_back, fit_amplitudes, model, cost)
end

function compute_lamp_fwhm_and_center_x(
    order::Int, λref::Float64, cfwhms::AbstractVector{T}, cxs::AbstractVector{T},
    λ::Float64, x) where {T<:Real}
    λpo = ((λ - λref) / λref) .^ (0:order)
    fwhm = λpo' * cfwhms
    center_x = λpo' * cxs
    sq_dist_to_center_x = (center_x - x)^2
    (fwhm, sq_dist_to_center_x)
end

encode_lamp_lkl_vmlmbvars(cfwhms::Vector{T}, cxs::Vector{T}) where {T<:Real} = hcat(cfwhms, cxs)


function decode_lamp_lkl_vmlmbvars(vmlmbvars::AbstractMatrix{<:Real})
    cfwhms = @view vmlmbvars[:, 1]
    cxs = @view vmlmbvars[:, 2]
    (cfwhms, cxs)
end

function (self::Lamp_LKL)(vmlmbvars::AbstractMatrix{<:Real})
    (cfwhms, cxs) = decode_lamp_lkl_vmlmbvars(vmlmbvars)
    (cost, _, _, _) = compute_lamp_cost_and_back_and_amplitudes(self, cfwhms, cxs)
    cost
end

function compute_lamp_cost_and_back_and_amplitudes(
    (; order, λref, bbox, data, lasers_pixels_λs)::Lamp_LKL,
    cfwhms::AbstractVector, cxs::AbstractVector
)

    lamp_image = compute_lamp_images(
        order, λref, cfwhms, cxs, lasers_pixels_λs, bbox)

    lamp_image_norm = lamp_image ./ sum(lamp_image; dims=1)

    (back, amplitudes...) = ChainRulesCore.@ignore_derivatives compute_lamp_backs_and_amplitudes(
        lamp_image_norm, data)

    model = @. (lamp_image_norm * amplitudes') + back

    cost = likelihood(data, model)
    #cost = likelihood(robustlikelihood(3.0), data, model)

    (cost, back, amplitudes, model)
end

compute_lamp_images((; λref, λs, fwhm_coefs, x_coefs, bbox)::LensletModel) =
    compute_lamp_images(length(fwhm_coefs) - 1, λref, fwhm_coefs, x_coefs, λs, bbox)

function compute_lamp_images(order::Int, λref::Float64, cfwhms::AbstractVector, cxs::AbstractVector,
    λs::Vector{<:AbstractFloat}, bbox::BoundingBox{Int}
)

    λpo = ((λs .- λref) ./ λref) .^ reshape(0:order, 1, order + 1)
    center_x = λpo * cxs

    fwhms = λpo * cfwhms

    (xs, ys) = axes(bbox)
    sq_dists = ((xs .- center_x') .^ 2)



    fwhm2sigma = 1 / (2 * sqrt(2 * log(2)))
    fw = -1 ./ (2 .* (fwhms .* fwhm2sigma) .^ 2)
    exp.(sq_dists .* reshape(fw, 1, :))
end
"""
    compute_lamp_backs_and_amplitudes(
        lamp_model::Matrix, data::Matrix, weights::Matrix) -> [background, amplitudes...]

Center and FWHM for each gaussian row is guessed by VMLMB. Optimal background and amplitudes can
be computed from them, this is what this function does.

# Arguments
- `lamp_model ::Matrix{Float64}`: of size `(W, H)`, each row is an image of a 2D gaussian
- `data ::Matrix`: of size `(W, H)`, contains lamp data for the lenslet
- `weights ::Matrix`: of size `(W, H)`, contains lamp data weights for the lenslet,
   high weight means high confidence, weight zero is for bad pixels.

# Returns
- `Vector{Float64}`: of size `(1+H)`, first value is the computed background, following values are
  the `H` computed amplitudes values, one for each row of the lamp model.

# Explanation

Caution: we call `W` and `H` the width and the height of the lenslet data, as we would see it
on a monitor screen. But since we store this in a Julia matrix, `W` is actually the number of rows
of the matrix, and `H` is actually the number of columns. This must be kept in mind while using
matrices operators.

if we define:
- `back` a scalar representing the background level in the lenslet
- `amps` a vector of size `height` containing amplitudes for each "gaussian row" of the lenslet
- `model` a matrix of size `(W, H)`, being `lamp_model` with each row multiplied by its
  respective amplitude, and additioned with the background:
  `model = (lamp_model .* ampsᵀ) .+ back`

the cost function (see `Lamp_LKL` is defined as:
`cost = sum(weights .* (model .- data).^2)`

if we define:
- `p` a vector of size `(1+H)`, containing `back`, followed by the `amps` values.
- `G` a matrix of size `(W*H, 1+H)`,
  we want `G⋅p` to be equal to `sum((lamp_model .* ampsᵀ) .+ back; dims=1)`,
  so: `G[1,:] .= 1`, `G[:,i] .= lamp_image[x,y]`
- `d` a vector of size `(W*H)` with `d[:] .= data[:,:]`
- `v` a vector of size `(W*H)` with `w[:] .= weights[:,:]`
- `V` a diagonal matrix of size `(W*H, W*H)` with `V[i,i] .= v[i]`
  
we can rewrite `cost` as:
`cost = (G⋅p - d)ᵀ ⋅ W ⋅ (G⋅amps - d)`

we want to find `amps` and `back` values that minimizes `cost`. So we derive `cost` by vector
`p`, and look at the expression when the derived `cost` equals zero.

first we rewrite `cost`:
`cost = (G⋅p)ᵀ⋅V⋅(G⋅p) + (dᵀ⋅V⋅d) - 2⋅dᵀ⋅V⋅G⋅p`

we derive by vector `p`:
`∂cost/∂p = (2⋅Gᵀ⋅V⋅G⋅p) - (2⋅dᵀ⋅V⋅G)` 

when this equals zero, we have an expression for `p`:
`(2⋅Gᵀ⋅V⋅G⋅p) - (2⋅dᵀ⋅V⋅G) = 0`
`p = (Gᵀ⋅V⋅G)⁻¹ ⋅ (dᵀ⋅V⋅G)`

so we define:
- `A = (Gᵀ⋅V⋅G)`, a matrix of size `(1+H, 1+H)`
- `b = (dᵀ⋅V⋅G)`, a vector of size `(1+H)`

which gives us:
`p = A⁻¹ ⋅ b`

In the function we compute `A⁻¹` and `b`, in a somehow efficient manner. A lots of values in `G`
are zeros so we avoid the basic matrix operation.
"""
function compute_lamp_backs_and_amplitudes(
    lamp_model::AbstractArray{T}, data::WeightedArray
) where {T<:Real}
    weights = get_precision(data)
    data = get_value(data)
    c = @. lamp_model * weights
    b = @. c * data
    a = @. lamp_model * c

    va = sum(a; dims=1)[:]
    vb = sum(b; dims=1)[:]
    vc = sum(c; dims=1)[:]

    za = (va .== T(0)) .|| (vb .<= T(0))

    va2 = map(i -> ifelse(za[i], T(1), va[i]), eachindex(va))
    vb2 = map(i -> ifelse(za[i], T(0), vb[i]), eachindex(vb))
    vc2 = map(i -> ifelse(za[i], T(0), vc[i]), eachindex(vc))

    # N = length(va2)
    # A ::Matrix{Float64}(undef,N+1,N+1)
    A = hcat(vcat(sum(weights), vc2), vcat(vc2', diagm(va2)))

    vb3 = vcat(sum(data .* weights), vb2[:])

    inv(A) * vb3
end

struct Profile{N}
    bbox::BoundingBox{Int64}
    ycenter::Float64
    cfwhm::Array{Float64,N}
    cx::Vector{Float64}
end

Profile(bbox::BoundingBox{Int}, cfwhm::AbstractArray, cx::AbstractVector) =
    Profile(bbox, mean(axes(bbox, 2)), cfwhm, cx)

((; bbox, ycenter, cfwhm, cx)::Profile)() = get_profile(bbox, ycenter, cfwhm, cx)


((; bbox, ycenter, cfwhm, cx)::Profile)(bbox2::BoundingBox{Int}) = get_profile(bbox2, ycenter, cfwhm, cx)

function get_profile(bbox::BoundingBox{Int64}, ycenter::Float64, cfwhm::Array{Float64,N}, cx::Vector{Float64}) where N


    xorder = length(cx)
    fwhmorder = size(cfwhm, 1)

    order = max(xorder, fwhmorder)

    ax, ay = axes(bbox)
    ypo = ((ay .- ycenter)) .^ reshape(0:order, 1, order + 1)

    xcenter = ypo[:, 1:xorder] * cx

    width = ypo[:, 1:fwhmorder] * cfwhm

    fwhm2sigma = 1 / (2 * sqrt(2 * log(2)))
    fw = @. -1 / (2 * (width * fwhm2sigma)^2)

    xc = (ax .- xcenter')
    if N == 1
        xc = (ax .- xcenter')
        dist = (xc .^ 2) .* reshape(fw, 1, :)
    elseif N == 2
        dist = min.(xc, 0) .^ 2 .* reshape(fw[:, 1], 1, :) .+ max.(xc, 0) .^ 2 .* reshape(fw[:, 2], 1, :)
    else
        error("get_profile : N must be 1 or 2")
    end


    img = exp.(dist)
    return img #./ sum(img; dims=1)
end



function extract_model(data::WeightedArray{T,N},
    profile::Profile;
    restrict=0.01,
    nonnegative=false,
    relative=false
) where {T,N}
    bbox = profile.bbox
    if relative
        (; value, precision) = data
    else
        if N > 2
            (; value, precision) = view(data, bbox, :)
        else
            (; value, precision) = view(data, bbox)
        end
    end
    model = profile()
    if restrict > 0
        model .*= (model .> restrict)
    end

    αprecision = dropdims(sum(model .^ 2 .* precision, dims=1), dims=1)
    α = dropdims(sum(model .* precision .* value, dims=1), dims=1) ./ αprecision

    nanpix = .!isnan.(α)
    if nonnegative
        positive = nanpix .& (α .>= T(0))
    else
        positive = nanpix
    end

    return WeightedArray(positive .* α, positive .* αprecision)
end

using OptimPackNextGen.Powell.Newuoa

function fit_profile(data::WeightedArray{T,N},
    profile::Profile{M};
    relative=false,
    optim=OptimParams()) where {T,N,M}

    fwhmorder = size(profile.cfwhm, 1)
    cxorder = length(profile.cx)
    if M == 1
        scale = vcat(10. .^ (-(1:(fwhmorder))), 10. .^ (-(1:(cxorder))))
    else
        scale = vcat(10. .^ (-(1:(fwhmorder))), 10. .^ (-(1:(cxorder))), 10. .^ (-(1:(cxorder))))
    end

    @unpack_OptimParams optim
    vec, re = Optimisers.destructure(profile)

    d = relative ? data : view(data, profile.bbox)
    f = build_loss(d, re)
    #f(x) = likelihood(ScaledL2Loss(dims=1, nonnegative=true), d, re(x)())
    #prep = prepare_gradient(f, ADbackend, vec)
    #fg!(x, grad) = DifferentiationInterface.value_and_gradient!(f, grad, prep, ADbackend, x)[1]
    #vmlmb!(fg!, vec; verb=verb, maxeval=maxeval, ftol=ftol, xtol=xtol, gtol=gtol, lower=lower, upper=upper)
    #    Newuoa.optimize!(f, vec, 1, 1e-3; scale=[1e-1, 1e-2, 1e-3, 1., 1e-2, 1e-2] .* ones(length(vec)), check=false, maxeval=10_000, verbose=0)
    Newuoa.optimize!(f, vec, 1, 1e-9; scale=scale, check=false, maxeval=10_000, verbose=0)
    #    @show f(vec)
    return re(vec)
end

build_loss(data, re) = x -> likelihood(ScaledL2Loss(dims=1, nonnegative=true), data, re(x)())

function get_meanx(data::WeightedArray{T,N}, bbox; relative=false) where {T,N}
    if relative
        (; value, precision) = data
    else
        (; value, precision) = view(data, bbox)
    end
    ax, ay = axes(bbox)

    return sum(value .* sqrt.(precision) .* ax) ./ sum(sqrt.(precision) .* value)
end

function refine_lamp_model(lamp, profiles, assigned_lenslets::AbstractVector{Bool}; loop=2, width=2)
    NLENS = length(profiles)
    lamp_profile = [WeightedArray(zeros(Float64, 40), zeros(Float64, 40)) for _ in 1:NLENS]
    model = zeros(Float64, size(lamp))
    tmodel = []
    progress = Progress(NLENS .* loop; showspeed=true)
    for _ ∈ 1:loop
        #res = lamp .- model
        res = WeightedArray(lamp.value .- model, lamp.precision)
        model = zeros(Float64, size(lamp))
        for i ∈ findall(assigned_lenslets)
            resi = WeightedArray(view(res, profiles[i].bbox).value .+ profiles[i]() .* reshape(lamp_profile[i].value, 1, :), view(res, profiles[i].bbox).precision)
            # resi = view(res,profiles[i].bbox) .+ profiles[i]() .* reshape(lamp_profile[i].value, 1, :)

            profiles[i] = fit_profile(resi, profiles[i]; relative=true)
            if any(isnan.(profiles[i].cfwhm))
                assigned_lenslets[i] = false
                continue
            end
            lamp_profile[i] = extract_model(resi, profiles[i]; relative=true)
            if any(isnan.(lamp_profile[i]))
                assigned_lenslets[i] = false
                continue
            end
            (; xmin, xmax, ymin, ymax) = profiles[i].bbox
            lbox = BoundingBox(xmin=xmin - width, xmax=xmax + width, ymin=ymin, ymax=ymax)
            p = profiles[i](lbox)
            view(model, lbox) .+= p .* reshape(lamp_profile[i].value, 1, :)
            next!(progress)
        end
        push!(tmodel, copy(model))
    end
    ProgressMeter.finish!(progress)
    return (; tmodel, lamp_profile, profiles)
end
using InterpolationKernels, SparseArrays

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



function find_index(knots::AbstractRange, sample)
    return (sample - first(knots)) / step(knots) + 1
end

function build_λrange(λs::AbstractMatrix{<:Real}; superres=1)
    nb_el = round(Int, size(λs, 1) * superres)
    blue = λs[1, :]
    red = λs[end, :]
    return range(start=minimum(blue), stop=maximum(red), step=median((red .- blue)) ./ nb_el)
end

function build_λrange(λs::Vector{Vector{Float64}}, assigned_lenslets; superres=1)
    idx = findall(assigned_lenslets)
    nb_el = round(Int, maximum(length.(λs[idx])) * superres)
    blue = minimum.(λs[idx])
    red = maximum.(λs[idx])
    return range(start=minimum(blue), stop=maximum(red), step=median((red .- blue)) ./ nb_el)
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


function buildAandB(knots, sampling, spectra, assigned_lenslets, α::T) where {T}
    nλ = length(knots)
    MI = Vector{SparseMatrixCSC{Float64,Int}}(undef, sum(assigned_lenslets))
    A = Matrix{Float64}(I, length(knots), length(knots))
    # b = Vector{Float64}(undef, length(knots))
    #A = zeros(Float64, nλ, nλ)
    b = zeros(Float64, nλ)

    for (i, idx) ∈ enumerate(findall(assigned_lenslets))
        MI[i] = build_sparse_interpolation_integration_matrix(knots, get_lower_uppersamples(sampling[idx])...)
        #   MI[i] = build_sparse_interpolation_matrix(knots, sampling[idx])
        if T <: Number
            (; value, precision) = ((1 ./ α) * spectra[idx])
        else
            (; value, precision) = ((1 ./ α[idx]) * spectra[idx])
        end

        b .+= Array(MI[i]' * (precision .* value))
        A .+= Array(MI[i]' * (precision .* MI[i]))
    end
    return MI, A, b
end

function spectral_refinement(coefs, data, spectrum_template, template_wavelength, reference_pixel)
    function loss(x)
        wvlngth = get_wavelength(x, reference_pixel, axes(data, 1))
        model = lamp_model(wvlngth, spectrum_template, template_wavelength)
        return likelihood(ScaledL2Loss(), data, model)
    end
    scale = 1e-6 .* vcat(10. .^ (-(1:length(coefs))))
    scale = 5.e-8 .* ones(length(coefs))

    coefs = copy(coefs)
    @show loss(coefs)
    Newuoa.optimize!(loss, coefs, 1, 1e-9; scale=scale, check=false, maxeval=10_000, verbose=1)
    @show loss(coefs)
    return coefs
end

function lamp_model(λ, spectrum_template, template_wavelength)
    lo, up = get_lower_uppersamples(λ)
    model = build_sparse_interpolation_integration_matrix(template_wavelength, lo, up) * spectrum_template
    return model
end



function laser_model(λ, fwhm_pixels, lasers_λs, data)
    idx = max.(2, [searchsortedlast(λ, l) for l ∈ lasers_λs])
    las = LaserModel(lasers_λs, fwhm_pixels .* (λ[idx] .- λ[idx.-1]))
    images = hcat(compute_laser_images(las, λ), ones(length(λ)))
    amplitude = compute_lasers_amplitudes(Val(length(lasers_λs) + 1), images, data)
    return images * amplitude
end

function spectral_refinement(coefs, lamp, lamp_template, wavelength, reference_pixel, lasers_λs, fwhm_pixels, laser)
    function loss(x)
        wvlngth = get_wavelength(x, reference_pixel, axes(lamp, 1))
        lamp_spectrum = lamp_model(wvlngth, lamp_template, wavelength)
        laser_spectrum = laser_model(wvlngth, fwhm_pixels, lasers_λs, laser)
        return likelihood(ScaledL2Loss(), lamp, lamp_spectrum) + likelihood(laser, laser_spectrum)
    end
    #scale = 1e-7 .* vcat(10. .^ (-(1:length(coefs))))
    scale = 1.e-8 .* ones(length(coefs))
    #coefs = copy(coefs)
    #   @show loss(coefs)
    Newuoa.optimize!(loss, coefs, 1, 1e-9; scale=scale, check=false, maxeval=10_000, verbose=0)
    #   @show loss(coefs)
    return coefs
end
using BandedMatrices

function estimate_template(λ, coefs, reference_pixel, spectra, assigned_lenslets)
    nλ = length(λ)
    transmission = zeros(Float64, length(assigned_lenslets))
    MI = Vector{SparseMatrixCSC{Float64,Int}}(undef, length(assigned_lenslets))
    A = zeros(Float64, nλ, nλ)
    diagA = 2 * ones(Float64, nλ)
    diagA[1] = 1
    diagA[end] = 1
    A = Array(BandedMatrix((0 => diagA, 1 => -1 * ones(nλ - 1), -1 => -1 * ones(nλ - 1)), (nλ, nλ)))

    b = zeros(Float64, nλ)
    foreach(findall(assigned_lenslets)) do idx
        (; value, precision) = spectra[idx]

        profile_wavelength = get_wavelength(coefs[idx], reference_pixel, 1:length(value))
        MI[idx] = build_sparse_interpolation_integration_matrix(λ, get_lower_uppersamples(profile_wavelength)...)
        b .+= Array(MI[idx]' * (precision .* value))
        A .+= Array(MI[idx]' * (precision .* MI[idx]))
    end


    template = A \ b

    OhMyThreads.tforeach(findall(assigned_lenslets)) do idx
        (; value, precision) = spectra[idx]
        m = (MI[idx] * template)
        transmission[idx] = sum((mp = m .* precision) .* value) / sum(m .* mp)
    end

    transmission .*= 1 ./ median(transmission[findall(assigned_lenslets)])

    return template, transmission
end

function recalibrate_wavelengths(λ,
    coefs,
    order,
    lamp_profile,
    laser_profile,
    lasers_λs,
    lasers_model,
    reference_pixel,
    assigned_lenslets;
    loop=2)

    template, transmission = estimate_template(λ, coefs, reference_pixel, lamp_profile, assigned_lenslets)

    new_coefs = similar(coefs)

    p = Progress(sum(assigned_lenslets) * loop; showspeed=true)

    for _ ∈ 1:loop
        @localize template @localize coefs OhMyThreads.tforeach(findall(assigned_lenslets)) do i
            if (order + 1) > length(coefs[i])
                coef = vcat(coefs[i], zeros(order - length(coefs[i]) + 1))
            else
                coef = copy(coefs[i])
            end
            try
                new_coefs[i] = spectral_refinement(coef, lamp_profile[i], template, λ, reference_pixel, lasers_λs, lasers_model[i].fwhm, laser_profile[i])
            catch e
                @warn "Spectral refinement failed for lenslet $i: $e"
                assigned_lenslets[i] = false
            end
            next!(p)
        end
        coefs = copy(new_coefs)

        template, transmission = estimate_template(λ, coefs, reference_pixel, lamp_profile, assigned_lenslets)

    end
    ProgressMeter.finish!(p)
    return coefs, template, transmission
end

function calibrate_profile(lamp,
    ; calib_params::PICParams=PICParams(),
    valid_lenslets::AbstractVector{Bool}=trues(calib_params.NLENS),
    loop=0,
    width=2
)


    @unpack_PICParams calib_params
    @unpack_BboxParams bbox_params

    size(valid_lenslets) == (NLENS,) || throw(ArgumentError("valid_lenslets must be of size NLENS"))


    bboxes = fill(BoundingBox{Int}(nothing), NLENS)
    profiles = Vector{Profile}(undef, NLENS)

    assigned_lenslets = falses(NLENS)

    @inbounds for i in findall(valid_lenslets)
        bbox = get_bbox(lasers_cxy0s_init[i, 1], lasers_cxy0s_init[i, 2]; bbox_params=bbox_params)
        if !ismissing(bbox)
            bboxes[i] = bbox
            assigned_lenslets[i] = true
            profiles[i] = Profile(bbox, lamp_cfwhms_init, vcat(PIC.get_meanx(lamp, bbox), zeros(profile_order)))
        end
    end

    profile_type = ZippedVector{WeightedValue{Float64},2,true,Tuple{Vector{Float64},Vector{Float64}}}
    lamp_profile = Vector{profile_type}(undef, NLENS)

    progress = Progress(sum(assigned_lenslets); showspeed=true)
    #Threads.@threads for i in findall(assigned_lenslets)
    # from https://discourse.julialang.org/t/optionally-multi-threaded-for-loop/81902/8?u=skleinbo
    _foreach = multi_thread ? OhMyThreads.tforeach : Base.foreach
    @allow_boxed_captures _foreach(findall(assigned_lenslets)) do i
        if sum(view(lamp, bboxes[i]).precision) == 0
            assigned_lenslets[i] = false
        else
            try

                profiles[i] = fit_profile(lamp, profiles[i])
                if any(isnan.(profiles[i].cfwhm))

                    throw("NaN found in profile for lenslet $i")
                end
                lamp_profile[i] = extract_model(lamp, profiles[i])

            catch e
                @debug "Error on lenslet $i" exception = (e, catch_backtrace())
                assigned_lenslets[i] = false
            end
        end
        next!(progress)
    end
    ProgressMeter.finish!(progress)

    model = zeros(Float64, size(lamp))

    progress = Progress(sum(assigned_lenslets) .* loop; showspeed=true)

    for _ ∈ 1:loop
        res = lamp .- model
        model = zeros(Float64, size(lamp))
        for i ∈ findall(assigned_lenslets)
            resi = WeightedArray(view(res, profiles[i].bbox).value .+ profiles[i]() .* reshape(lamp_profile[i].value, 1, :), view(res, profiles[i].bbox).precision)
            # resi = view(res,profiles[i].bbox) .+ profiles[i]() .* reshape(lamp_profile[i].value, 1, :)

            profiles[i] = fit_profile(resi, profiles[i]; relative=true)
            if any(isnan.(profiles[i].cfwhm))
                assigned_lenslets[i] = false
                continue
            end
            lamp_profile[i] = extract_model(resi, profiles[i]; relative=true)
            if any(isnan.(lamp_profile[i]))
                assigned_lenslets[i] = false
                continue
            end
            (; xmin, xmax, ymin, ymax) = profiles[i].bbox
            lbox = BoundingBox(xmin=xmin - width, xmax=xmax + width, ymin=ymin, ymax=ymax)
            p = profiles[i](lbox)
            view(model, lbox) .+= p .* reshape(lamp_profile[i].value, 1, :)
            next!(progress)
        end
    end
    ProgressMeter.finish!(progress)

    return profiles, bboxes, assigned_lenslets, lamp_profile, model
end

function spectral_calibration(
    lasers,
    lamp_profiles,
    profiles;
    assigned_lenslets=trues(length(profiles)),
    calib_params::PICParams=PICParams(),
    loop=2,
    superres=1,
    final_spectral_order=3
)

    @unpack_PICParams calib_params
    @unpack_BboxParams bbox_params

    profile_type = ZippedVector{WeightedValue{Float64},2,true,Tuple{Vector{Float64},Vector{Float64}}}
    laser_profile = Vector{profile_type}(undef, NLENS)
    laser_model = LaserModel([7.0, 20.0, 35.0], [2.0, 2.0, 2.0])
    coefs = Vector{Vector{Float64}}(undef, NLENS)
    λ = Vector{Vector{Float64}}(undef, NLENS)
    las = Vector{typeof(laser_model)}(undef, NLENS)

    #Threads.@threads for i in findall(assigned_lenslets)
    # from https://discourse.julialang.org/t/optionally-multi-threaded-for-loop/81902/8?u=skleinbo
    _foreach = multi_thread ? (@localize coefs OhMyThreads.tforeach) : Base.foreach
    progress = Progress(sum(assigned_lenslets); showspeed=true)
    @localize coefs _foreach(findall(assigned_lenslets)) do i
        if sum(view(lasers, profiles[i].bbox).precision) == 0
            assigned_lenslets[i] = false
        else
            try
                laser_profile[i] = extract_model(lasers, profiles[i])
                las[i] = fit_laser(laser_profile[i], laser_model)

                if std(las[i].position .- laser_model.position) > 1
                    throw("Laser position too far from initial guess for lenslet $i")
                end
                W = get_laser_precision(las[i], laser_profile[i])
                if any(diag(W) .< 1e-4)
                    throw("W singular  for lenslet $i")
                end
                coefs[i] = spectral_calibration(spectral_order, reference_pixel, lasers_λs, las[i].position, W)
                if any(isnan.(coefs[i]))
                    throw("NaN found in coefs for lenslet $i")
                end
                λ[i] = get_wavelength(coefs[i], reference_pixel, axes(laser_profile[i], 1))

            catch e
                @debug "Error on lenslet $i" exception = e
                assigned_lenslets[i] = false
            end
        end
        next!(progress)
    end
    finish!(progress)
    lλ = build_λrange(λ, assigned_lenslets; superres=superres)

    coefs, template, transmission = recalibrate_wavelengths(
        lλ,
        coefs,
        final_spectral_order,
        lamp_profiles,
        laser_profile,
        lasers_λs,
        las,
        reference_pixel,
        assigned_lenslets;
        loop=loop)
    return coefs, template, transmission, lλ, las, laser_profile, assigned_lenslets
end
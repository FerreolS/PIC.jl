
"""
    Lasers_LKL(model::LensletModel,wavelengths::AbstractArray{<:Real,1},data::AbstractArray,weight::AbstractArray)

Build the likelihood function for a given lenslet
* `lmodel`: model of the lenslet
* `laser`: wavelengths of the illumination lasers
* `data` : data
* `weight`: precision (ie inverse variance) of the data
"""
@concrete struct Lasers_LKL
    nλ::Int
    order::Int  # order of the polynomial
    lasers_λs <: AbstractVector
    λref <: Real   # reference wavelength
    bbox::BoundingBox{Int}
    data <: WeightedArray

end

function encode_lasers_lkl_vmlmbvars(
    fwhms::Vector{T}, cxs::Vector{T}, cys::Vector{T}
) where {T<:AbstractFloat}
    nλ = length(fwhms)
    length(cxs) == length(cys) || throw(ArgumentError)
    length(cxs) == nλ || throw(ArgumentError)
    vmlmbvars = similar(fwhms, 3 * nλ)
    vmlmbvars[1:nλ] = fwhms
    vmlmbvars[(nλ+1):2nλ] = cxs
    vmlmbvars[(2nλ+1):3nλ] = cys
    return vmlmbvars
end

function decode_lasers_lkl_vmlmbvars(
    nλ::Int, vmlmbvars::AbstractVector
)
    fwhms = vmlmbvars[1:nλ]
    cxs = vmlmbvars[(nλ+1):2nλ]
    cys = vmlmbvars[(2nλ+1):3nλ]
    return (fwhms, cxs, cys)
end

function (self::Lasers_LKL)(vmlmbvars::AbstractVector)
    (fwhms, cxs, cys) = decode_lasers_lkl_vmlmbvars(self.nλ, vmlmbvars)
    (cost, _, _) = compute_lasers_cost_and_amplitudes(self, cxs, cys, fwhms)
    cost
end

function fit_lens_lasers(
    lasers_lkl::Lasers_LKL,
    lasers_fwhms_init::AbstractVector,
    lasers_cxs_init::AbstractVector,
    lasers_cys_init::AbstractVector;
    optim=OptimParams()
)
    @unpack_OptimParams optim
    vmlmbvars = encode_lasers_lkl_vmlmbvars(lasers_fwhms_init, lasers_cxs_init, lasers_cys_init)


    if true
        #vmlmb!(lasers_lkl, vmlmbvars; verb=false, ftol=(0.0, 1e-8), maxeval=500, autodiff=true)
        grad = similar(vmlmbvars)
        prep = prepare_gradient(lasers_lkl, ADbackend, vmlmbvars)
        fg!(x, grad) = DifferentiationInterface.value_and_gradient!(lasers_lkl, grad, prep, ADbackend, x)[1]
        vmlmb!(fg!, vmlmbvars; verb=verb, maxeval=maxeval, ftol=ftol, xtol=xtol, gtol=gtol, lower=lower, upper=upper)
        #xopt, info = prima(lasers_lkl, vmlmbvars; maxfun=10_000, ftarget=length(lens_lasers_data))

    else
        f(x) = lasers_lkl(x)
        OptimPackNextGen.Powell.Newuoa.newuoa!(f, vmlmbvars, 1.e-1, 1.e-6, verbose=false, maxeval=5000)
    end
    (fit_fwhms, fit_cxs, fit_cys) = decode_lasers_lkl_vmlmbvars(lasers_lkl.nλ, vmlmbvars)

    (cost, fit_amplitudes, model) = compute_lasers_cost_and_amplitudes(
        lasers_lkl, fit_cxs, fit_cys, fit_fwhms)

    (fit_cxs, fit_cys, fit_fwhms, fit_amplitudes, model, cost)
end


function compute_laser_center(
    order::Int, λref::Float64, cxs::AbstractVector, cys::AbstractVector, λ::Float64
)
    #λpo = ((λ - λref) / λref) .^ (1:order)
    #center_x = cxs[1] + sum(cxs[2:end] .* λpo)
    #center_y = cys[1] + sum(cys[2:end] .* λpo)
    λpo = ((λ - λref) / λref) .^ (0:order)
    center_x = cxs' * λpo
    center_y = cys' * λpo
    (center_x, center_y)
end


"""
    GaussianModel2(fwhm::Float64, x::AbstractArray)

Compute the value at lenslets_coords sqrt(r) 1D centered Gaussian
* `fwhm` : full-width at half maximum
* `x`:  squared sampled lenslets_coords

Equivalent to `GaussianModel(1.,fwhm, sqrt(x))`
"""
function GaussianModel2(fwhm::Real, x::T) where {T<:Real}
    fwhm2sigma = 1 / (2 * sqrt(2 * log(2)))
    exp(-x / (2 * (fwhm * fwhm2sigma)^2))
end

function GaussianModel2(t::NTuple{2,T}) where {T<:Real}
    return GaussianModel2(t[1], t[2])
end


function compute_lasers_images(
    ::Val{N}, order::Int, λref::Float64, cxs::AbstractVector{T}, cys::AbstractVector{T},
    fwhms::AbstractVector{T}, lasers_λs::Vector{<:AbstractFloat}, bbox::BoundingBox{Int}
) where {N,T}

    λpo = ((lasers_λs .- λref) ./ λref) .^ reshape(0:order, 1, order + 1)
    center_x = λpo * cxs
    center_y = λpo * cys


    (xs, ys) = axes(bbox)
    a = ((xs .- center_x') .^ 2)
    b = ((ys .- center_y') .^ 2)


    sq_dists = reshape(a, :, 1, N) .+ reshape(b, 1, :, N)

    fwhm2sigma = T(1 / (2 * sqrt(2 * log(2))))
    fw = -T(1) ./ (T(2) .* (fwhms .* fwhm2sigma) .^ 2)
    exp.(sq_dists .* reshape(fw, 1, 1, N))
end

function compute_lasers_cost_and_amplitudes(
    (; nλ, order, λref, bbox, lasers_λs, data)::Lasers_LKL,
    cxs::Vector{T}, cys::Vector{T}, fwhms::Vector{T}) where {T<:Real}

    laser_images = compute_lasers_images(Val(nλ), order, λref, cxs, cys, fwhms, lasers_λs, bbox)

    amplitudes = ChainRulesCore.@ignore_derivatives compute_lasers_amplitudes(Val(nλ), laser_images, data)

    #model = sum(i -> laser_images[i] .* amplitudes[i], 1:lkl.nλ)
    #model = sum(laser_images .* amplitudes)
    model = reshape(reshape(laser_images, :, nλ) * amplitudes, size(bbox))
    #model = mapreduce(x -> (.*)(x...), +, zip(laser_images, amplitudes))
    #model =amplitudes' * laser_images

    cost = likelihood(data, model)
    #cost = likelihood(data, model, likelihoodfunc=robustlikelihood(3.0))
    (cost, amplitudes, model)
end

"""
    compute_lasers_amplitudes(
        lasers_models::Vector{Matrix}, data::Matrix, weights::Matrix) -> amplitudes::Vector

Center and FWHM of each laser spot is guessed by VMLMB. Optimal amplitude can be computed from
them, this is what this function does.

# Arguments
- `lasers_models ::Vector{Matrix{Float64}}`: containing `nλ` matrices of size `(W, H)`, each containing
  a gaussian model for a laser spot, without background and with amplitude `1`
- `data ::Matrix`: of size `(W, H)`, containing lasers data, for the lenslet bbox
- `weights ::Matrix`: of size `(W, H)`, containing lasers data weights for the lenslet
  bbox, high weight means high confidence, weight zero is for bad pixels

# Returns
- `Vector{Float64}`: of size `nλ`, the computed amplitude for each laser

# Explanation

Caution: we call `W` and `H` the width and the height of the lenslet data, as we would see it
on a monitor screen. But since we store this in a Julia matrix, `W` is actually the number of rows
of the matrix, and `H` is actually the number of columns. This must be kept in mind while using
matrices operators.

if we define:
- `amps` a vector of size `(nλ)` containing amplitudes for each gaussian laser spot
- `sum_models` a matrix of size `(W, H)`, the sum of lasers models multiplied by their
  respective amplitude:
  `sum_models = sum(lasers_models .* amps; dims=3)`

the cost function (see `Lasers_LKL`) is defined as:
`cost = sum(weights .* (sum_models .- data).^2)`

if we define:
- `G` a matrix of size `(W*H, nλ)` with `G[:,i] .= lasers_models[:,:,i]`
- `d` a vector of size `(W*H)` with `d[:] .= data[:,:]`
- `v` a vector of size `(W*H)` with `w[:] .= weights[:,:]`
- `V` a diagonal matrix of size `(W*H, W*H)` with `V[i,i] .= v[i]`

we can rewrite `cost` as:
`cost = (G⋅amps - d)ᵀ ⋅ V ⋅ (G⋅amps - d)`

we want to find the `amps` value that minimizes `cost`. So we derive `cost` by vector `amps`,
and look at the expression when the derived `cost` equals zero.

First we rewrite `cost`:
`cost = (G⋅amps)ᵀ⋅V⋅(G⋅amps) + (dᵀ⋅V⋅d) - 2⋅dᵀ⋅V⋅G⋅amps`

we derive by vector `amp`:
`∂cost/∂amps = (2⋅Gᵀ⋅V⋅G⋅amps) - (2⋅dᵀ⋅V⋅G)`

when this equals zero, we have an expression for `amps`:
`(2⋅Gᵀ⋅V⋅G⋅amps) - (2⋅dᵀ⋅V⋅G) = 0`
`amps = (Gᵀ⋅V⋅G)⁻¹ ⋅ (dᵀ⋅V⋅G)`

so we define:
- `A = (Gᵀ⋅V⋅G)`, a matrix of size `(nλ, nλ)`
- `b = (dᵀ⋅V⋅G)`, a vector of size `(nλ)`

which gives us:
`amps = A⁻¹ ⋅ b`

In the function we compute `A⁻¹` and `b`.
"""
function compute_lasers_amplitudes(::Val{N},
    lasers_models::Vector{Matrix{T}}, d::WeightedArray
) where {T<:Real,N}
    data = get_value(d)
    precision = WeightedData.get_precision(d)
    model = [reshape(lasers_models[i], :) for i in 1:N]
    d = view(data, :)
    w = view(precision, :)

    A = @MMatrix zeros(T, N, N)
    b = @MVector zeros(T, N)


    @inbounds for index = 1:N
        mw = model[index] .* w
        b[index] = mw' * d
        A[index, index] = mw' * model[index]
        for i = 1:index-1
            A[i, index] = A[index, i] = mw' * model[i]
        end
    end
    return inv(A) * b

end

function compute_lasers_amplitudes(::Val{N},
    lasers_models::Array{T,3}, d::WeightedArray
) where {T<:Real,N}
    data = get_value(d)
    precision = WeightedData.get_precision(d)
    d = view(data, :)
    w = view(precision, :)


    if false
        model = reshape(lasers_models, :, N)

        mw = (w .* model)

        A = mw' * model
        b = mw' * d
    else
        model = [reshape(lasers_models[:, :, i], :) for i in 1:N]
        A = @MMatrix zeros(T, N, N)
        b = @MVector zeros(T, N)

        @inbounds for index = 1:N
            mw = model[index] .* w
            b[index] = mw' * d
            A[index, index] = mw' * model[index]
            for i = 1:index-1
                A[i, index] = A[index, i] = mw' * model[i]
            end
        end
    end
    return inv(A) * b

end

function compute_lasers_dists_and_λmap!(
    λrange::AbstractVector{Float64}, bbox::BoundingBox{Int}, lasers_order::Int, λref::Float64,
    laser_cxs::AbstractVector, laser_cys::AbstractVector,
    laser_pixels_dists::AbstractMatrix, laser_pixels_λs::AbstractMatrix
)

    previous_index = 0
    I0 = first(CartesianIndices(bbox)) - CartesianIndex(1, 1)
    for I in CartesianIndices(bbox)
        previous_index = max(1, previous_index - 5)
        for (index, λ) in enumerate(λrange[previous_index:end])
            (laser_cx, laser_cy) = compute_laser_center(
                lasers_order, λref, laser_cxs, laser_cys, λ)
            dist_to_laser_x = (I[1] - laser_cx)
            dist_to_laser = sqrt(dist_to_laser_x^2 + (I[2] - laser_cy)^2)
            r = sign(dist_to_laser_x) * dist_to_laser
            if isnan(laser_pixels_dists[I-I0]) || abs(r) < abs(laser_pixels_dists[I-I0])
                laser_pixels_dists[I-I0] = r
                laser_pixels_λs[I-I0] = λ
            else
                break
            end
            previous_index += 1
        end
    end

end


function compute_lasers_λmap!(
    λrange::AbstractVector{Float64}, bbox::BoundingBox{Int}, order::Int, λref::Float64,
    laser_cxs::AbstractVector, laser_cys::AbstractVector,
    laser_pixels_dists::AbstractVector, laser_pixels_λs::AbstractVector
)
    x = round(Int, middle(axes(bbox, 1)))
    previous_index = 0

    for (ny, y) in enumerate(axes(bbox, 2))
        previous_index = max(1, previous_index - 5)
        for (index, λ) in enumerate(λrange[previous_index:end])

            λpo = ((λ - λref) / λref) .^ (0:order)
            laser_cy = laser_cys' * λpo
            dist_to_laser = abs(y - laser_cy)
            if isnan(laser_pixels_dists[ny]) || abs(dist_to_laser) < abs(laser_pixels_dists[ny])
                laser_pixels_dists[ny] = dist_to_laser
                laser_pixels_λs[ny] = λ
            else
                break
            end
            previous_index += 1
        end
    end

end

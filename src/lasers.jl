
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
    (cost, amplitudes) = compute_lasers_cost_and_amplitudes(self, cxs, cys, fwhms)
    cost
end

function fit_lens_lasers(
    lasers_lkl::Lasers_LKL,
    lasers_fwhms_init::AbstractVector,
    lasers_cxs_init::AbstractVector,
    lasers_cys_init::AbstractVector
)

    vmlmbvars = encode_lasers_lkl_vmlmbvars(lasers_fwhms_init, lasers_cxs_init, lasers_cys_init)


    #vmlmb!(lasers_lkl, vmlmbvars; verb=false, ftol=(0.0, 1e-8), maxeval=500, autodiff=true)
    grad = similar(vmlmbvars)
    prep = prepare_gradient(lasers_lkl, AutoZygote(), vmlmbvars)
    fg!(x, grad) = DifferentiationInterface.value_and_gradient!(lasers_lkl, grad, prep, AutoZygote(), x)[1]
    vmlmb!(fg!, vmlmbvars; verb=false, ftol=(0.0, 1e-8), maxeval=500)
    #xopt, info = prima(lasers_lkl, vmlmbvars; maxfun=10_000, ftarget=length(lens_lasers_data))
    (fit_fwhms, fit_cxs, fit_cys) = decode_lasers_lkl_vmlmbvars(lasers_lkl.nλ, vmlmbvars)

    (cost, fit_amplitudes) = compute_lasers_cost_and_amplitudes(
        lasers_lkl, fit_cxs, fit_cys, fit_fwhms)

    (fit_cxs, fit_cys, fit_fwhms, fit_amplitudes)
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

function compute_laser_image(
    laser_center_x, laser_center_y, fwhm, bbox::BoundingBox{Int}
)
    (xs, ys) = axes(bbox)
    sq_dists = ((xs .- laser_center_x) .^ 2) .+ ((ys .- laser_center_y) .^ 2)'
    GaussianModel2.(fwhm, sq_dists)
end

@generated function compute_lasers_images(
    ::Val{N}, order::Int, λref::Float64, cxs::AbstractVector, cys::AbstractVector,
    fwhms::AbstractVector, lasers_λs::Vector{<:AbstractFloat}, bbox::BoundingBox{Int}
) where {N}
    if N != 3 && N != 4
        throw(ArgumentError("compute_lasers_images: N must be 3 or 4"))
    end
    code = Expr(:block)
    for i in 1:N
        push!(code.args, quote
            λpo = ((lasers_λs[$i] - λref) / λref) .^ (0:order)
            $(Symbol("laser_center_x" * "$i")) = λpo' * cxs
            $(Symbol("laser_center_y" * "$i")) = λpo' * cys
            #  laser_center_y$i = λpo' * cys
        end)
    end
    #(laser_center_x, laser_center_y) = compute_laser_center(order, λref, cxs, cys, lasers_λs[i])
    if N == 3
        push!(code.args, quote
            return [compute_laser_image(laser_center_x1, laser_center_y1, fwhms[1], bbox),
                compute_laser_image(laser_center_x2, laser_center_y2, fwhms[2], bbox),
                compute_laser_image(laser_center_x3, laser_center_y3, fwhms[3], bbox)]
        end)
    elseif N == 4
        push!(code.args, quote
            return [compute_laser_image(laser_center_x1, laser_center_y1, fwhms[1], bbox),
                compute_laser_image(laser_center_x2, laser_center_y2, fwhms[2], bbox),
                compute_laser_image(laser_center_x3, laser_center_y3, fwhms[3], bbox),
                compute_laser_image(laser_center_x4, laser_center_y4, fwhms[4], bbox)]
        end)
    end
    return code
end

function compute_lasers_cost_and_amplitudes(
    lkl::Lasers_LKL, cxs::Vector{T}, cys::Vector{T}, fwhms::Vector{T}
) where {T<:Real}

    laser_images = compute_lasers_images(
        Val(lkl.nλ), lkl.order, lkl.λref, cxs, cys, fwhms, lkl.lasers_λs, lkl.bbox)

    amplitudes = ChainRulesCore.@ignore_derivatives compute_lasers_amplitudes(Val(lkl.nλ), laser_images, lkl.data)

    #model = sum(i -> laser_images[i] .* amplitudes[i], 1:lkl.nλ)
    model = sum(laser_images .* amplitudes)
    #model = mapreduce(x -> (.*)(x...), +, zip(laser_images, amplitudes))
    #model =amplitudes' * laser_images

    cost = likelihood(lkl.data, model)

    (cost, amplitudes)
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
    data = get_data(d)
    precision = get_precision(d)
    model = [reshape(lasers_models[i], :) for i in 1:N]
    d = view(data, :)
    w = view(precision, :)

    A = @MMatrix zeros(T, N, N)
    b = @MVector zeros(T, N)

    mw = similar(model)

    @inbounds for index = 1:N
        mw[index] = model[index] .* w
        b[index] = mw[index]' * d
        A[index, index] = mw[index]' * model[index]
        for i = 1:index-1
            A[i, index] = A[index, i] = mw[index]' * model[i]
        end
    end
    return inv(A) * b

end

function compute_lasers_dists_and_λmap!(
    λrange::AbstractVector{Float64}, bbox::BoundingBox{Int}, lasers_order::Int, λref::Float64,
    laser_cxs::Vector{<:AbstractFloat}, laser_cys::Vector{<:AbstractFloat},
    laser_pixels_dists::AbstractMatrix{<:AbstractFloat}, laser_pixels_λs::AbstractMatrix{<:AbstractFloat}
)::Nothing

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

    nothing
end

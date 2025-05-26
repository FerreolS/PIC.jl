
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
    data = get_data(data)
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

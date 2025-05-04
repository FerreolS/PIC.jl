function fit_lens_lamp(
    bbox ::BoundingBox{Int},
    lens_lamp_data        ::AbstractMatrix{<:Real},
    lens_lamp_weights     ::AbstractMatrix{<:Real},
    lens_lasers_pixels_λs ::AbstractMatrix{<:Real},
    ; lamp_order ::Int,
      λref ::Float64,
      lamp_cfwhms_init ::Vector{Float64},
      lamp_cxs_init    ::Vector{Float64}
) ::Tuple{Vector{Float64},Vector{Float64},Float64,Vector{Float64}}

    lamp_lkl = Lamp_LKL(lamp_order, λref, bbox, lens_lasers_pixels_λs,
                                lens_lamp_data, lens_lamp_weights)

    vmlmbvars = encode_lamp_lkl_vmlmbvars(lamp_cfwhms_init, lamp_cxs_init)

    vmlmb!(lamp_lkl, vmlmbvars; verb=false, ftol=(0.0,1e-8), maxeval=500, autodiff=true)
    
    (fit_cfwhms, fit_cxs) = decode_lamp_lkl_vmlmbvars(vmlmbvars)

    (cost, fit_back, fit_amplitudes) = compute_lamp_cost_and_back_and_amplitudes(
        lamp_lkl, fit_cfwhms, fit_cxs)

    (fit_cfwhms, fit_cxs, fit_back, fit_amplitudes)
end

function compute_lamp_fwhm_and_center_x(
    order::Int, λref::Float64, cfwhms::AbstractVector{Float64}, cxs::AbstractVector{Float64},
    λ::Float64, x::Int
) ::NTuple{2,Float64}
    λpo = ((λ - λref) / λref).^(1:order)
    fwhm = cfwhms[1] + sum(cfwhms[2:end] .* λpo)
    center_x = cxs[1] + sum(cxs[2:end] .* λpo)
    sq_dist_to_center_x = (center_x - x)^2
    (fwhm, sq_dist_to_center_x)
end

struct Lamp_LKL{L<:AbstractMatrix{<:Real},D<:AbstractMatrix{<:Real},W<:AbstractMatrix{<:Real}}
    order::Int
    λref::Float64
    bbox::BoundingBox{Int}
    lasers_pixels_λs::L
    data::D
    weights::W
    function Lamp_LKL{L,D,W}(order, λref, bbox, lasers_pixels_λs, data, weights) where {L,D,W}
        size(bbox) == size(lasers_pixels_λs) || throw(ArgumentError)
        size(bbox) == size(data)             || throw(ArgumentError)
        size(bbox) == size(weights)          || throw(ArgumentError)
        new{L,D,W}(order, λref, bbox, lasers_pixels_λs, data, weights)
    end
end

function Lamp_LKL(
    order::Int, λref::Float64, bbox::BoundingBox{Int}, lasers_pixels_λs::L, data::D, weights::W
) where {L<:AbstractMatrix{<:Real},D<:AbstractMatrix{<:Real},W<:AbstractMatrix{<:Real}}
    Lamp_LKL{L,D,W}(order, λref, bbox, lasers_pixels_λs, data, weights)
end

function encode_lamp_lkl_vmlmbvars(cfwhms::Vector{Float64}, cxs::Vector{Float64}) ::Matrix{Float64}
    vmlmbvars = hcat(cfwhms, cxs)
end

function decode_lamp_lkl_vmlmbvars(vmlmbvars::Matrix{Float64}) ::NTuple{2,AbstractVector{Float64}}
    cfwhms = view(vmlmbvars,:,1)
    cxs    = view(vmlmbvars,:,2)
    (cfwhms, cxs)
end

function (self::Lamp_LKL)(vmlmbvars::Matrix{Float64}) ::Float64
    (cfwhms, cxs) = decode_lamp_lkl_vmlmbvars(vmlmbvars)
    (cost, back, amps) = compute_lamp_cost_and_back_and_amplitudes(self, cfwhms, cxs)
    cost
end

function compute_lamp_cost_and_back_and_amplitudes(
    lkl::Lamp_LKL, cfwhms::AbstractVector{Float64}, cxs::AbstractVector{Float64}
) ::Tuple{Float64,Float64,Vector{Float64}}
    
    bbox_rx = axes(lkl.bbox, 1)
    
    lamp_image = [
        GaussianModel2(compute_lamp_fwhm_and_center_x(
            lkl.order, lkl.λref, cfwhms, cxs, lkl.lasers_pixels_λs[x,y], bbox_rx[x]))
        for x in 1:size(lkl.bbox,1), y in 1:size(lkl.bbox,2) ]

    lamp_image_norm = lamp_image ./ sum(lamp_image; dims=1)

    (back, amplitudes...) = compute_lamp_backs_and_amplitudes(
        lamp_image_norm, lkl.data, lkl.weights)

    model = @. (lamp_image_norm * amplitudes') + back

    cost = sum(@. lkl.weights * (lkl.data - model)^2 )
    
    (cost, back, amplitudes)
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
    lamp_model, data::D, weights::W
) ::Vector{Float64} where {D<:AbstractMatrix{<:Real},
                                                   W<:AbstractMatrix{<:Real}}
    
    c = @. lamp_model * weights
    b = @. lamp_model * data * weights
    a = @. lamp_model^2 * weights
    
    va = sum(a; dims=1)[:]
    vb = sum(b; dims=1)[:]
    vc = sum(c; dims=1)[:]
    
    za = (va .== 0) .|| (vb .<= 0)

    va2 = map(i -> za[i] ? 1 : va[i], eachindex(va))
    vb2 = map(i -> za[i] ? 0 : vb[i], eachindex(vb))
    vc2 = map(i -> za[i] ? 0 : vc[i], eachindex(vc))

    # N = length(va2)
    # A ::Matrix{Float64}(undef,N+1,N+1)
    A = hcat( vcat(sum(weights), vc2), vcat(vc2', diagm(va2)) )

    vb3 = vcat(sum(data .* weights), vb2[:])

    inv(A) * vb3
end

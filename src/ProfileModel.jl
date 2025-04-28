mutable struct ProfileModel
    λref::Float64   # reference wavelength
    order::Int  # order of the polynomial
    cλs::Vector{Float64} # coefficients of the polynomial along the wavelength axis
    cxs::Vector{Float64} # coefficients of the polynomial along the x axis
    function ProfileModel(λref,order,cλs,cxs)
        order ≥ 0               || throw(ArgumentError)
        length(cλs) == (order+1) || throw(ArgumentError)
        length(cxs) == (order+1) || throw(ArgumentError)
        new(λref, order, cλs, cxs)
    end
end

function ProfileModel(λref::Float64, order::Int)
    cλs = zeros(order+1)
    cxs = zeros(order+1)
    cλs[1] = 1
    cxs[1] = 1
    ProfileModel(λref, order, cλs, cxs)
end

function ProfileModel(λref::Float64, coefs::Vector{Float64})
    order = Int(length(coefs) / 2) - 1
    cλs = coefs[1 : (order+1)]
    cxs = coefs[(order+2) : end]
    ProfileModel(λref, order, cλs, cxs)
end

function compute_lamp_fwhm_and_center_x(
    order::Int, λref::Float64, cλs::AbstractVector{Float64}, cxs::AbstractVector{Float64},
    λ::Float64, x::Int
) ::NTuple{2,Float64}
    λpo = ((λ - λref) / λref).^(1:order)
    fwhm = cλs[1] + sum(cλs[2:end] .* λpo)
    center_x = cxs[1] + sum(cxs[2:end] .* λpo)
    sq_dist_to_center_x = (center_x - x)^2
    (fwhm, sq_dist_to_center_x)
end

function (self::ProfileModel)(λ::Float64, x::Int) ::NTuple{2,Float64} # [?, pix]
    λpo = ((λ-self.λref)/self.λref).^(1:self.order)
    w = self.cλs[1] + sum(self.cλs[2:end] .* λpo)
    gaussian_cxs = self.cxs[1] + sum(self.cxs[2:end] .* λpo) # [pix coord]
    dist_to_gaussian_cxs = (gaussian_cxs - x)^2             # [pix]
    return (w, dist_to_gaussian_cxs)
end

function updateProfileModel!(
    self::ProfileModel, cλs::AbstractVector{Float64}, cxs::AbstractVector{Float64}
) ::Nothing
    length(cλs) == (self.order+1) || throw(ArgumentError)
    length(cxs) == (self.order+1) || throw(ArgumentError)
    self.cλs = cλs
    self.cxs = cxs
    nothing
end

struct Lamp_LKL{D<:AbstractMatrix{<:Real},W<:AbstractMatrix{<:Real}}
    order::Int
    λref::Float64
    bbox::BoundingBox{Int}
    lasers_pixels_λs::AbstractMatrix{Float64}
    data::D
    weights::W
    function Lamp_LKL{D,W}(order, λref, bbox, lasers_pixels_λs, data, weights) where {D,W}
        size(bbox) == size(lasers_pixels_λs) || throw(ArgumentError)
        size(bbox) == size(data)             || throw(ArgumentError)
        size(bbox) == size(weights)          || throw(ArgumentError)
        new{D,W}(order, λref, bbox, lasers_pixels_λs, data, weights)
    end
end

function Lamp_LKL(
    order::Int, λref::Float64, bbox::BoundingBox{Int},
    lasers_pixels_λs::AbstractMatrix{Float64}, data::D, weights::W
) where {D<:AbstractMatrix{<:Real},W<:AbstractMatrix{<:Real}}
    Lamp_LKL{D,W}(order, λref, bbox, lasers_pixels_λs, data, weights)
end

function encode_lamp_lkl_fitvars(cλs::Vector{Float64}, cxs::Vector{Float64}) ::Matrix{Float64}
    fitvars = [ cλs ;; cxs ]
end

function decode_lamp_lkl_fitvars(fitvars::Matrix{Float64}) ::NTuple{2,AbstractVector{Float64}}
    cλs = view(fitvars,:,1)
    cxs = view(fitvars,:,2)
    (cλs, cxs)
end

function (self::Lamp_LKL)(fitvars::Matrix{Float64}) ::Float64
    (cλs, cxs) = decode_lamp_lkl_fitvars(fitvars)
    (cost, back, amps) = compute_lamp_back_and_amps(self, cλs, cxs)
    cost
end

function compute_lamp_back_and_amps(
    lkl::Lamp_LKL, cλs::AbstractVector{Float64}, cxs::AbstractVector{Float64}
) ::Tuple{Float64,Float64,Vector{Float64}}
    
    bbox_rx = axes(lkl.bbox, 1)
    
    lamp_image = [
        GaussianModel2(compute_lamp_fwhm_and_center_x(
            lkl.order, lkl.λref, cλs, cxs, lkl.lasers_pixels_λs[x,y], bbox_rx[x]))
        for x in 1:size(lkl.bbox,1), y in 1:size(lkl.bbox,2) ]

    lamp_image_norm = lamp_image ./ sum(lamp_image; dims=1)

    (back, amps) = compute_lamp_backs_and_amps(lamp_image_norm, lkl.data, lkl.weights)

    model = @. (lamp_image_norm * amps') + back

    cost = sum(@. lkl.weights * (lkl.data - model)^2 )
    
    (cost, back, amps)
end

"""
    compute_lamp_backs_and_amps(
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
function compute_lamp_backs_and_amps(
    lamp_model, data::D, weights::W
) ::Tuple{Float64,Vector{Float64}} where {D<:AbstractMatrix{<:Real},
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

    N = length(va2)
    A = Matrix{Float64}(undef,N+1,N+1)
    A = hcat( vcat(sum(weights), vc2), vcat(vc2', diagm(va2)) )

    vb3 = vcat(sum(data .* weights), vb2[:])

    (back, amplitudes...) = inv(A) * vb3
end

"""
    Lasers_LKL(model::LensletModel,wavelengths::AbstractArray{<:Real,1},data::AbstractArray,weight::AbstractArray)

Build the likelihood function for a given lenslet
* `lmodel`: model of the lenslet
* `laser`: wavelengths of the illumination lasers
* `data` : data
* `weight`: precision (ie inverse variance) of the data
"""
struct Lasers_LKL{D<:AbstractMatrix{<:Real},W<:AbstractMatrix{<:Real}}
    nλ::Int
    order::Int64  # order of the polynomial
    lasers_λs::Vector{Float64}
    λref::Float64   # reference wavelength
    bbox::BoundingBox{Int}
    data::D
    weights::W
    function Lasers_LKL{D,W}(
        nλ, order, lasers_λs, λref, bbox, data, weights
    ) where {D,W}
        length(lasers_λs) == nλ        || throw(ArgumentError)
        size(data)       == size(bbox) || throw(ArgumentError)
        size(weights)    == size(bbox) || throw(ArgumentError)
        new{D,W}(nλ, order, lasers_λs, λref, bbox, data, weights)
    end
end

function Lasers_LKL(
    nλ::Int, order::Int64, lasers_λs::Vector{Float64}, λref::Float64, bbox::BoundingBox{Int},
    data::D, weights::W
) where {D<:AbstractMatrix{<:Real},W<:AbstractMatrix{<:Real}}
    Lasers_LKL{D,W}(nλ, order, lasers_λs, λref, bbox, data, weights)
end

function encode_lasers_lkl_fitvars(
    fwhms::Vector{Float64}, cxs::Vector{Float64}, cys::Vector{Float64}
) ::Vector{Float64}
    fitvars = Float64[]
    append!(fitvars, fwhms)
    for i in 1:length(cxs)
        push!(fitvars, cxs[i])
        push!(fitvars, cys[i])
    end
    fitvars
end

function decode_lasers_lkl_fitvars(
    nλ::Int, fitvars::Vector{Float64}
) ::NTuple{3,Vector{Float64}}
    fwhms = fitvars[1:nλ]
    cxs  = fitvars[ (nλ+1) : 2 : (end-1) ]
    cys  = fitvars[ (nλ+2) : 2 :  end    ]
    (fwhms, cxs, cys)
end

function (self::Lasers_LKL)(fitvars::Vector{Float64}) ::Float64
    (fwhms, cxs, cys) = decode_lasers_lkl_fitvars(self.nλ, fitvars)
    (cost, amplitudes) = compute_lasers_cost_and_amplitudes(self, cxs, cys, fwhms)
    cost
end

function compute_laser_center(
    order::Int, λref::Float64, cxs::Vector{Float64}, cys::Vector{Float64}, λ::Float64
) ::NTuple{2,Float64}
    λpo = ((λ - λref)/λref).^(1:order)
    x = cxs[1] + sum(cxs[2:end] .* λpo)
    y = cys[1] + sum(cys[2:end] .* λpo)
    (x, y)
end

function compute_laser_image(
    laser_center_x::Float64, laser_center_y::Float64, fwhm::Float64, bbox::BoundingBox{Int}
) ::Matrix{Float64}
    (xs, ys) = axes(bbox)
    sq_dists = ((xs .- laser_center_x).^2) .+ ((ys .- laser_center_y).^2)'
    GaussianModel2.(fwhm, sq_dists)
end

function compute_lasers_images(
    nλ::Int, order::Int, λref::Float64, cxs::Vector{Float64}, cys::Vector{Float64},
    fwhms::Vector{Float64}, lasers_λs::Vector{Float64}, bbox::BoundingBox{Int}
) ::Vector{Matrix{Float64}}
    map(1:nλ) do i
        (laser_center_x, laser_center_y) = compute_laser_center(order, λref, cxs, cys, lasers_λs[i])
        matrix = compute_laser_image(laser_center_x, laser_center_y, fwhms[i], bbox)
    end
end

function compute_lasers_cost_and_amplitudes(
    lkl::Lasers_LKL, cxs::Vector{Float64}, cys::Vector{Float64}, fwhms::Vector{Float64}
) ::Tuple{Float64,Vector{Float64}}

    laser_images = compute_lasers_images(
        lkl.nλ, lkl.order, lkl.λref, cxs, cys, fwhms, lkl.lasers_λs, lkl.bbox)

    amplitudes = compute_lasers_amplitudes(laser_images, lkl.data, lkl.weights)
    
    model = sum(i -> laser_images[i] .* amplitudes[i], 1:lkl.nλ)
    
    cost = sum(lkl.weights .* (lkl.data .- model).^2)
    
    (cost, amplitudes)
end

 """
    compute_amplitude(images::Array{Float64,3}, data::Matrix, weights::Matrix) -> Vector{Float64}

From a gaussian laser images model, and data and weights from the lasers file, compute the
amplitude for each gaussian laser spot, for a lenslet.

# Arguments
- `images` is an `Array{Float64,3}` of size `(W,H,nλ)`, containing the gaussian model for each
  laser spot, each without background and with theoretical integral equal to `1`.
- `data` is a matrix of size `(W,H)` containing lasers data, for the lenslet bbox
- `weights` is a matrix of size `(W,H)` containing lasers data weights, for the lenslet bbox.
  high weight means high confidence, weight zero is for bad pixels

if we define:
- `(W,H)` as the size of the bbox of the lenslet
- `nλ` as the number of laser images
- `amp` as a vector of size `nλ` containing amplitudes for each gaussian laser spot
- `model` as the sum of images multiplied by their respective amplitude:
  `model = sum(images .* amp; dims=3)`

the cost function (see `Lasers_LKL`) is defined as:
`cost = sum(weights .* (model .- data).^2)`

if we define:
- `G` a matrix of size `(W*H, nλ)` with `G[:,i] .= images[:,:,i]`
- `d` a vector of size `(W*H)` with `d[:] .= data[:,:]`
- `w` a vector of size `(W*H)` with `w[:] .= weights[:,:]`
- `W` a diagonal matrix of size `(W*H,W*H)` with `W[i,i] .= w[i]`

we can rewrite `cost` as:
`cost = (G⋅amp - d)ᵀ ⋅ W ⋅ (G⋅amp - d)`

we want to find the `amp` value that minimizes `cost`. So we derive `cost` by vector `amp`,
and look at the expression when the derived `cost` equals zero.

First we rewrite `cost`:
`cost = (G⋅amp)ᵀ⋅W⋅(G⋅amp) + (dᵀ⋅W⋅d) - 2⋅dᵀ⋅W⋅G⋅amp`

we derive by vector `amp`:
`∂cost/∂amp = (2⋅Gᵀ⋅W⋅G⋅amp) - (2⋅dᵀ⋅W⋅G)`

when this equals zero, we have an expression for `amp`:
`(2⋅Gᵀ⋅W⋅G⋅amp) - (2⋅dᵀ⋅W⋅G) = 0`
`amp = (Gᵀ⋅W⋅G)⁻¹ ⋅ (dᵀ⋅W⋅G)`

finally we define:
- `A = (Gᵀ⋅W⋅G)`, a matrix of size `(nλ,nλ)`
- `b = (dᵀ⋅W⋅G)`, a vector of size `(nλ)`
We compute `A` and `b` in the function, inverse `A`, then we have a value for `amp`.
"""
function compute_lasers_amplitudes(
    images::Vector{Matrix{Float64}}, data::AbstractMatrix, weights::AbstractMatrix
) ::Vector{Float64}
    
    A = [ sum(images[i] .* weights .* images[j]) for i in 1:3, j in 1:3 ]
    b = [ sum(data .* weights .* images[i]) for i in 1:3 ]
    
    amp = inv(A) * b
end


function compute_lasers_dists_and_λmap!(
    λrange::AbstractVector{Float64}, bbox::BoundingBox{Int}, lasers_order::Int, λref::Float64,
    laser_cxs::Vector{Float64}, laser_cys::Vector{Float64},
    laser_pixels_dists::AbstractMatrix{Float64}, laser_pixels_λs::AbstractMatrix{Float64}
) ::Nothing

    previous_index = 0
    I0 = first(CartesianIndices(bbox)) - CartesianIndex(1,1)
    for I in CartesianIndices(bbox)
        previous_index = max(1, previous_index-5)
        for (index,λ) in enumerate(λrange[previous_index:end])
            (laser_cx, laser_cy) = compute_laser_center(
                lasers_order, λref, laser_cxs, laser_cys, λ)
            dist_to_laser_x = (I[1] - laser_cx)
            dist_to_laser = sqrt(dist_to_laser_x^2 + (I[2] - laser_cy)^2)
            r = sign(dist_to_laser_x) * dist_to_laser
            if isnan(laser_pixels_dists[I-I0]) || abs(r) < abs(laser_pixels_dists[I-I0])
                laser_pixels_dists[I-I0] = r;
                laser_pixels_λs[I-I0] = λ;
            else
                break
            end
            previous_index += 1
        end
    end
    
    nothing
end
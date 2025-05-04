function fit_lens_lasers(
    bbox ::BoundingBox{Int},
    lens_lasers_data    ::AbstractMatrix{<:Real},
    lens_lasers_weights ::AbstractMatrix{<:Real},
    ; nλ ::Int,
      lasers_λs ::Vector{Float64},
      λref ::Float64,
      lasers_order ::Int,
      lasers_fwhms_init ::Vector{Float64},
      lasers_cxs_init   ::Vector{Float64},
      lasers_cys_init   ::Vector{Float64}
) ::NTuple{4,Vector{Float64}}
      
    lasers_lkl = Lasers_LKL(
        nλ, lasers_order, lasers_λs, λref, bbox, lens_lasers_data, lens_lasers_weights)
    
    vmlmbvars = encode_lasers_lkl_vmlmbvars(lasers_fwhms_init, lasers_cxs_init, lasers_cys_init)

    vmlmb!(lasers_lkl, vmlmbvars; verb=false, ftol=(0.0,1e-8), maxeval=500, autodiff=true)
    
    (fit_fwhms, fit_cxs, fit_cys) = decode_lasers_lkl_vmlmbvars(nλ, vmlmbvars)
        
    (cost, fit_amplitudes) = compute_lasers_cost_and_amplitudes(
        lasers_lkl, fit_cxs, fit_cys, fit_fwhms)

    (fit_cxs, fit_cys, fit_fwhms, fit_amplitudes)
end

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

function encode_lasers_lkl_vmlmbvars(
    fwhms::Vector{Float64}, cxs::Vector{Float64}, cys::Vector{Float64}
) ::Vector{Float64}
    vmlmbvars = Float64[]
    append!(vmlmbvars, fwhms)
    for i in 1:length(cxs)
        push!(vmlmbvars, cxs[i])
        push!(vmlmbvars, cys[i])
    end
    vmlmbvars
end

function decode_lasers_lkl_vmlmbvars(
    nλ::Int, vmlmbvars::Vector{Float64}
) ::NTuple{3,Vector{Float64}}
    fwhms = vmlmbvars[1:nλ]
    cxs  = vmlmbvars[ (nλ+1) : 2 : (end-1) ]
    cys  = vmlmbvars[ (nλ+2) : 2 :  end    ]
    (fwhms, cxs, cys)
end

function (self::Lasers_LKL)(vmlmbvars::Vector{Float64}) ::Float64
    (fwhms, cxs, cys) = decode_lasers_lkl_vmlmbvars(self.nλ, vmlmbvars)
    (cost, amplitudes) = compute_lasers_cost_and_amplitudes(self, cxs, cys, fwhms)
    cost
end

function compute_laser_center(
    order::Int, λref::Float64, cxs::Vector{Float64}, cys::Vector{Float64}, λ::Float64
) ::NTuple{2,Float64}
    λpo = ((λ - λref) / λref).^(1:order)
    center_x = cxs[1] + sum(cxs[2:end] .* λpo)
    center_y = cys[1] + sum(cys[2:end] .* λpo)
    (center_x, center_y)
end

"""
    GaussianModel2(fwhm::Float64, x::AbstractArray)

Compute the value at lenslets_coords sqrt(r) 1D centered Gaussian
* `fwhm` : full-width at half maximum
* `x`:  squared sampled lenslets_coords

Equivalent to `GaussianModel(1.,fwhm, sqrt(x))`
"""
function GaussianModel2(fwhm::T, x::T) ::T where {T<:Real}
    fwhm2sigma = 1 / (2 * sqrt(2 * log(2)))
    exp(-x / (2 * (fwhm * fwhm2sigma)^2))
end

function GaussianModel2(t::NTuple{2,T}) ::T where {T<:Real} return GaussianModel2(t[1], t[2]) end

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
    
    cost = sum(@. lkl.weights * (lkl.data - model)^2)
    
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
function compute_lasers_amplitudes(
    lasers_models::Vector{Matrix{Float64}}, data::AbstractMatrix, weights::AbstractMatrix
) ::Vector{Float64}
    
    A = [ sum(lasers_models[i] .* weights .* lasers_models[j]) for i in 1:3, j in 1:3 ]
    b = [ sum(data .* weights .* lasers_models[i]) for i in 1:3 ]
    
    amps = inv(A) * b
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



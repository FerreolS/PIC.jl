"""
    LasersModel(λ0::Float64,order::Int32,cxs::Array{Float64,1},cys::Array{Float64,1})

The lasers model giving the position of a wavelength on the detector
* `λ0` is the reference wavelength
* `order` is the order of the polynomials
* `cxs` is an array of coefficients of the polynomial along the x axis
* `cys` is an array of coefficients of the polynomial along the y axis
"""
struct LasersModel
    nλ::Int
    order::Int64  # order of the polynomial
    λref::Float64   # reference wavelength
    cxs::Vector{Float64} # coefficients of the polynomial along the x axis
    cys::Vector{Float64} # coefficients of the polynomial along the y axis
    fwhms::Vector{Float64}
    amplitudes::Vector{Float64}
    function LasersModel(nλ, order, λref, cxs, cys, fwhms, amplitudes)
        nλ ≥ 2                   || throw(ArgumentError)
        order ≥ 1                || throw(ArgumentError)
        length(cxs) == (order+1) || throw(ArgumentError)
        length(cys) == (order+1) || throw(ArgumentError)
        length(fwhms) == nλ      || throw(ArgumentError)
        length(amplitudes) == nλ || throw(ArgumentError)
        new(nλ, order, λref, cxs, cys, fwhms, amplitudes)
    end
end

function LasersModel(nλ::Int, order::Int, λref::Float64)
    nλ ≥ 2    || throw(ArgumentError)
    order ≥ 1 || throw(ArgumentError)
    cxs = Vector{Float64}(undef, order+1)
    cys = Vector{Float64}(undef, order+1)
    fwhms = Vector{Float64}(undef, nλ)
    amplitudes = Vector{Float64}(undef, nλ)
    LasersModel(nλ, order, λref, cxs, cys, fwhms, amplitudes)
end

function compute_λ_peak(
    order::Int, λref::Float64, cxs::Vector{Float64}, cys::Vector{Float64}, λ::Float64
) ::NTuple{2,Float64}
    λpo = ((λ - λref)/λref).^(1:order)
    x = cxs[1] + sum(cxs[2:end] .* λpo)
    y = cys[1] + sum(cys[2:end] .* λpo)
    (x, y)
end

function compute_λ_peak(lasers_model::LasersModel, λ::Float64) ::NTuple{2,Float64}
    compute_λ_peak(lasers_model.order, lasers_model.λref, lasers_model.cxs, lasers_model.cys, λ)
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
    bbox::BoundingBox{Int}
    lasers_model::LasersModel
    lasers_λs::Vector{Float64}
    data::D
    weights::W
    spots::Array{Float64,3}
    function Lasers_LKL{D,W}(
        nλ, bbox, lasers_model, lasers_λs, data, weights, spots
    ) where {D,W}
        length(lasers_λs) == nλ        || throw(ArgumentError)
        nλ > lasers_model.order          || throw(ArgumentError)
        size(data)       == size(bbox) || throw(ArgumentError)
        size(weights)    == size(bbox) || throw(ArgumentError)
        size(spots)[1:2] == size(bbox) || throw(ArgumentError)
        size(spots,3) == nλ            || throw(ArgumentError)
        new{D,W}(nλ, bbox, lasers_model, lasers_λs, data, weights, spots)
    end
end

function Lasers_LKL(
    bbox::BoundingBox{Int}, lasers_model::LasersModel, lasers_λs::Vector{Float64},
    data::D, weights::W
) where {D<:AbstractMatrix{<:Real},W<:AbstractMatrix{<:Real}}
    nλ = length(lasers_λs)
    spots = zeros(Float64, size(bbox)..., nλ)
    Lasers_LKL{D,W}(nλ, bbox, lasers_model, lasers_λs, data, weights, spots)
end

function encode_lasers_lkl_fitvars(
    fwhm::Vector{Float64}, cxs::Vector{Float64}, cys::Vector{Float64}
) ::Vector{Float64}
    length(cxs) == length(cys) || throw(ArgumentError)
    fitvars = Float64[]
    append!(fitvars, fwhm)
    for i in 1:length(cxs)
        push!(fitvars, cxs[i])
        push!(fitvars, cys[i])
    end
    fitvars
end

function decode_lasers_lkl_fitvars(
    nλ::Int, order::Int, fitvars::Vector{Float64}
) ::NTuple{3,Vector{Float64}}
    length(fitvars) == (nλ + 2 * (order + 1)) || throw(ArgumentError)
    fwhm = fitvars[1:nλ]
    cxs  = fitvars[ (nλ+1) : 2 : (end-1) ]
    cys  = fitvars[ (nλ+2) : 2 :  end    ]
    (fwhm, cxs, cys)
end

function (self::Lasers_LKL)(fitvars::Vector{Float64}) ::Float64

    (fwhm, cxs, cys) = decode_lasers_lkl_fitvars(self.nλ, self.lasers_model.order, fitvars)

    (xs,ys) = axes(self.bbox) # extracting bounding box range
    
    spots_buffer = Zygote.Buffer(self.spots)
    @inbounds for (index,λ) in enumerate(self.lasers_λs)  # For all laser
#        (mx, my) = compute_λ_peak(self.lasers_model, λ)  # center of the index-th Gaussian spot
        (mx, my) = compute_λ_peak(self.lasers_model.order, self.lasers_model.λref, cxs, cys, λ)
        xys = ((xs .- mx).^2) .+ ((ys .- my).^2)'
        spots_buffer[:,:,index] = GaussianModel2.(fwhm[index], xys)
    end
    spots = copy(spots_buffer)
    amplitudes = compute_amplitudes(spots, self.data, self.weights)
    sumspot = zeros(Float64, size(self.bbox))
    @inbounds for i in 1:self.nλ
        if isnan(amplitudes[i])
            @debug "NaN amplitude"
        else
            sumspot += amplitudes[i] * spots[:,:,i]
        end
    end
    
    Zygote.@ignore begin
        self.lasers_model.cxs .= cxs
        self.lasers_model.cys .= cys
        self.lasers_model.fwhms .= fwhm
        self.lasers_model.amplitudes .= amplitudes
    end
    
    return sum(self.weights .* (self.data .- sumspot).^2)
end

 """
    compute_amplitude(spots::Array{Float64,3}, data::Matrix, weights::Matrix) -> Vector{Float64}

From a gaussian laser spots model, and data and weights from the lasers file, compute the
amplitude for each gaussian laser spot, for a lenslet.

# Arguments
- `spots` is an `Array{Float64,3}` of size `(W,H,nλ)`, containing the gaussian model for each
  laser spot, each without background and with theoretical integral equal to `1`.
- `data` is a matrix of size `(W,H)` containing lasers data, for the lenslet bbox
- `weights` is a matrix of size `(W,H)` containing lasers data weights, for the lenslet bbox.
  high weight means high confidence, weight zero is for bad pixels

if we define:
- `(W,H)` as the size of the bbox of the lenslet
- `nλ` as the number of laser spots
- `amp` as a vector of size `nλ` containing amplitudes for each gaussian laser spot
- `model` as the sum of spots multiplied by their respective amplitude:
  `model = sum(spots .* amp; dims=3)`

the cost function (see `Lasers_LKL`) is defined as:
`cost = sum(weights .* (model .- data).^2)`

if we define:
- `G` a matrix of size `(W*H, nλ)` with `G[:,i] .= spots[:,:,i]`
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
function compute_amplitudes(
    spots::AbstractArray{Float64,3}, data::AbstractMatrix, weights::AbstractMatrix
) ::Vector{Float64}
    
    A = [ sum(spots[:,:,i] .* weights .* spots[:,:,j]) for i in 1:3, j in 1:3 ]
    b = [ sum(data .* weights .* spots[:,:,i]) for i in 1:3 ]
    
    amp = inv(A) * b
end

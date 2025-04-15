"""
    DispersionModel(λ0::Float64,order::Int32,cx::Array{Float64,1},cy::Array{Float64,1})

The dispersion model giving the position of a wavelength on the detector
* `λ0` is the reference wavelength
* `order` is the order of the polynomials
* `cx` is an array of coefficients of the polynomial along the x axis
* `cy` is an array of coefficients of the polynomial along the y axis
"""
mutable struct DispersionModel
    λref::Float64   # reference wavelength
    order::Int64  # order of the polynomial
    cx::Vector{Float64} # coefficients of the polynomial along the x axis
    cy::Vector{Float64} # coefficients of the polynomial along the y axis
    function DispersionModel(λref, order, cx, cy)
        order ≥ 0               || throw(ArgumentError)
        length(cx) == (order+1) || throw(ArgumentError)
        length(cy) == (order+1) || throw(ArgumentError)
        new(λref, order, cx, cy)
    end
end

function DispersionModel(λref::Float64, order::Int)
    cx = zeros(order+1)
    cx[1]=1
    cy = zeros(order+1)
    cy[1]=1
    DispersionModel(λref,order,cx,cy)
end

"""
    (self::DispersionModel)(λ::Float64)

compute the position `(x,y)`  of the wavelength `λ`
according to the dispersion law `DispersionModel`.

### Example
```
D = DispersionModel(λ0, order, cx, cy);
(x,y) = D(λ)
```
"""
function (self::DispersionModel)(λ::Float64)
    λpo = ((λ - self.λref)/self.λref).^(1:self.order)
    x = self.cx[1] + sum(self.cx[2:end] .* λpo)
    y = self.cy[1] + sum(self.cy[2:end] .* λpo)
    (x, y)
end


"""
    updateDispersionModel!(self::DispersionModel, cxs::Vector{Float64}, cys::Vector{Float64}) -> Nothing

Update the coefficients  of the DispersionModel .
* `self`: DispersionModel object
* `cxs` : vector containing the X polynomial coefficients
* `cys` : vector containing the X polynomial coefficients
"""
function updateDispersionModel!(self::DispersionModel, cxs::Vector{Float64}, cys::Vector{Float64}) ::Nothing
    length(cxs) == (self.order+1) || throw(ArgumentError)
    length(cys) == (self.order+1) || throw(ArgumentError)
    self.cx = cxs
    self.cy = cys
    nothing
end

"""
    Dispersion_LKL(model::LensletModel,wavelengths::AbstractArray{<:Real,1},data::AbstractArray,weight::AbstractArray)

Build the likelihood function for a given lenslet
* `lmodel`: model of the lenslet
* `laser`: wavelengths of the illumination lasers
* `data` : data
* `weight`: precision (ie inverse variance) of the data
"""
struct Dispersion_LKL{D<:AbstractMatrix{<:Real},W<:AbstractMatrix{<:Real}}
    nλ::Int
    bbox::BoundingBox{Int}
    disp_model::DispersionModel
    lasers_λs::Vector{Float64}
    data::D
    weights::W
    spots::Array{Float64,3}
    amplitude::Vector{Float64}
    function Dispersion_LKL{D,W}(
        nλ, bbox, disp_model, lasers_λs, data, weights, spots, amplitude
    ) where {D,W}
        length(lasers_λs) == nλ        || throw(ArgumentError)
        size(data) == size(weights)    || throw(ArgumentError)
        size(spots,3) == nλ            || throw(ArgumentError)
        length(amplitude) == nλ        || throw(ArgumentError)
        nλ > disp_model.order          || throw(ArgumentError)
        size(spots)[1:2] == size(bbox) || throw(ArgumentError)
        new{D,W}(nλ, bbox, disp_model, lasers_λs, data, weights, spots, amplitude)
    end
end

function Dispersion_LKL(
    bbox::BoundingBox{Int}, disp_model::DispersionModel, lasers_λs::Vector{Float64}, data::D, weights::W
) where {D<:AbstractMatrix{<:Real},W<:AbstractMatrix{<:Real}}
    nλ = length(lasers_λs)
    spots = zeros(Float64, size(bbox)..., nλ)
    amplitude = zeros(Float64, nλ)
    Dispersion_LKL{D,W}(nλ, bbox, disp_model, lasers_λs, data, weights, spots, amplitude)
end

function (self::Dispersion_LKL)(xs::Vector{Float64}) ::Float64

    fwhm = xs[1:self.nλ]
    cxs = xs[ (self.nλ+1) : 2 : (end-1) ]
    cys = xs[ (self.nλ+2) : 2 :  end    ]
    updateDispersionModel!(self.disp_model, cxs, cys)
    
    (xs,ys) = axes(self.bbox) # extracting bounding box range
    
    spots_buffer = Zygote.Buffer(self.spots)
    @inbounds for (index,λ) in enumerate(self.lasers_λs)  # For all laser
        (mx, my) = self.disp_model(λ)  # center of the index-th Gaussian spot
        xys = ((xs .- mx).^2) .+ ((ys .- my).^2)'
        spots_buffer[:,:,index] = GaussianModel2.(fwhm[index], xys)
    end
    spots = copy(spots_buffer)
    Zygote.@ignore self.amplitude .= compute_amplitude(spots, self.data, self.weights)
    sumspot = zeros(Float64, size(self.bbox))
    @inbounds for i in 1:self.nλ
        sumspot += self.amplitude[i] * spots[:,:,i]
    end
    return sum(self.weights .* (self.data .- sumspot).^2)
end

 """
    compute_amplitude(spots::Array{Float64,3}, data::Matrix, weights::Matrix) -> Vector{Float64}

From a gaussian laser spots model, and data and weights from the dispersion file, compute the
amplitude for each gaussian laser spot, for a lenslet.

# Arguments
- `spots` is an `Array{Float64,3}` of size `(W,H,nλ`), containing the gaussian model for each
  laser spot, each without background and with theoretical integral equal to `1`.
- `data` is a matrix of size `(W,H)` containing dispersion data, for the lenslet bbox
- `weights` is a matrix of size `(W,H)` containing dispersion data weights, for the lenslet bbox.
  high weight means high confidence, weight zero is for bad pixels

if we define:
- `(W,H)` as the size of the bbox of the lenslet
- `nλ` as the number of laser spots
- `amp` as a vector of size `nλ` containing amplitudes for each gaussian laser spot
- `model` as the sum of spots multiplied by their respective amplitude:
  `model = sum(spots .* amp; dims=3)`

the cost function (see `Dispersion_LKL`) is defined as:
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
function compute_amplitude(
    spots::AbstractArray{Float64,3}, data::AbstractMatrix, weights::AbstractMatrix
) ::Vector{Float64}

    size(data) == size(weights) == size(spots)[1:2] || throw(ArgumentError)
    
    nλ = size(spots,3)
    A = @MMatrix zeros(Float64,nλ,nλ)
    b = @MVector zeros(Float64,nλ)

    @inbounds for i in 1:nλ
        for j in 1:i
            A[i,j] = A[j,i] = sum(spots[:,:,i] .* weights .* spots[:,:,j])
        end
        b[i] = sum(data .* weights .* spots[:,:,i])
    end
    
    amp = inv(A) * b
end

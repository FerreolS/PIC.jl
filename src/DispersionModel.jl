"""
    DispersionModel(λ0::Float64,order::Int32,cx::Array{Float64,1},cy::Array{Float64,1})

The dispersion model giving the position of a wavelength on the detector
* `λ0` is the reference wavelength
* `order` is the order of the polynomials
* `cx` is an array of coefficients of the polynomial along the x axis
* `cy` is an array of coefficients of the polynomial along the y axis
"""
mutable struct DispersionModel
    λ0::Float64   # reference wavelength
    order::Int64  # order of the polynomial
    cx::Vector{Float64} # coefficients of the polynomial along the x axis
    cy::Vector{Float64} # coefficients of the polynomial along the y axis
    function DispersionModel(λ0, order, cx, cy)
        order ≥ 0               || throw(ArgumentError)
        length(cx) == (order+1) || throw(ArgumentError)
        length(cy) == (order+1) || throw(ArgumentError)
        new(λ0, order, cx, cy)
    end
end

function DispersionModel(λ0::Float64, order::Int)
    cx = zeros(order+1)
    cx[1]=1
    cy = zeros(order+1)
    cy[1]=1
    DispersionModel(λ0,order,cx,cy)
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
    λpo = ((λ - self.λ0)/self.λ0).^(1:self.order)
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
    Zygote.@ignore self.amplitude .= updateAmplitude(self.nλ, spots, self.data, self.weights)
    sumspot = zeros(Float64, size(self.bbox))
    @inbounds for i in 1:self.nλ
        sumspot += self.amplitude[i] * spots[:,:,i]
    end
    return sum(self.weights .* (self.data .- sumspot).^2)
end

 """
        updateAmplitude(nλ,m,d,W)

    return the `nλ` amplitudes `a` according the the model `m`, the data and the precision `W`
    such that
    `a = argmin_a || a*m - D||^2_W`
    where
    * `nλ` : is the number of spots in the model
    * `m`:  is the model composed of `nλ` images of spots
    * `d`:  is the data
    * `W`: is the precision (inverse variance) of the data
"""
function updateAmplitude(
    nλ::Int, spots::AbstractArray{Float64,3}, data::D, weights::W
) ::Vector{Float64} where {D<:AbstractMatrix{<:Real},W<:AbstractMatrix{<:Real}}
    A = @MMatrix zeros(Float64,nλ,nλ)
    b = @MVector zeros(Float64,nλ)
    spots_buffer = similar(spots)
    @inbounds for i in 1:nλ
        spots_buffer[:,:,i] .= spots[:,:,i] .* weights
        b[i] = sum(spots_buffer[:,:,i] .* data)
        A[i,i] = sum(spots_buffer[:,:,i] .* spots[:,:,i])
        for j in 1:(i-1)
            A[j,i] = A[i,j] = sum(spots_buffer[:,:,i] .* spots[:,:,j])
        end
    end
    return inv(A) * b
end

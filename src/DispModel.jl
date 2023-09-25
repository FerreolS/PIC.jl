
"""
    DispModel(λ0::Float64,order::Int32,cx::Array{Float64,1},cy::Array{Float64,1})

The dispersion model giving the position of a wavelength on the detector
* `λ0` is the reference wavelength
* `order` is the order of the polynomials
* `cx` is an array of coefficients of the polynomial along the x axis
* `cy` is an array of coefficients of the polynomial along the y axis
"""
mutable struct DispModel
    λ0::Float64   # reference wavelength
    order::Int64  # order of the polynomial
    cx::Array{Float64,1} # coefficients of the polynomial along the x axis
    cy::Array{Float64,1} # coefficients of the polynomial along the y axis

    function DispModel(λ0,order::Int64,cλ::Vector{Float64},cy::Vector{Float64})
        @assert length(cλ)==(order+1) "coefficients size does not match the order"
        @assert length(cy)==(order+1) "coefficients size does not match the order"
        new(λ0,order,cλ,cy)
    end
end



function DispModel(λ0,order::Int64)
    cx = zeros(order+1)
    cx[1]=1
    cy = zeros(order+1)
    cy[1]=1
    DispModel(λ0,order,cx,cy)
end

"""
    (self::DispModel)(λ::Float64)

compute the position `(x,y)`  of the wavelength `λ`
according to the dispersion law `DispModel`.

### Example
```
D = DispModel(λ0, order, cx, cy);
(x,y) = D(λ)
```
"""
function (self::DispModel)(λ::Float64)
    # x = self.cx[1];
    # y = self.cy[1];
    # for o in 1:self.order
    #     λpo = (( λ - self.λ0)/self.λ0 )^(o)
    #     x += self.cx[o + 1]  * λpo;
    #     y += self.cy[o + 1]  * λpo;
    # end

    λpo = (( λ - self.λ0)/self.λ0 ).^(1:self.order)
    x = self.cx[1] +sum(self.cx[2:end] .* λpo)
    y = self.cy[1] +sum(self.cy[2:end] .* λpo)
    
    return (x, y)
end

"""
    UpdateDispModel(self::DispModel, C::Array{Float64,2})

Update the coefficients  of the DispModel .
* `self`: DispModel object
* `C` : array containing the polynomial coefficients.
"""
function UpdateDispModel(self::DispModel, C::Array{Float64,2})
    @assert size(C)==(2,self.order+1) "coefficients array does not have the right size"
    self.cx = C[1,:];
    self.cy = C[2,:];
    return self
end

"""
    LaserModel(nλ::Int,λlaser::Array{Float64,1},amplitude::Array{Float64,1},fwhm::Array{Float64,1})

Model of the laser illumination.
It consist on:
* `nλ` the number of laser
* `λlaser` an array of the wavelengths of the lasers
* `amplitude` an array of the amplitude of the maximum of the Gaussian spot
* `fwhm` an array of the full width at half maximum of the Gaussians
"""
mutable struct LaserModel
    nλ::Int64
    λlaser::Array{Float64,1}# wavelength of the laser
    amplitude::Array{Float64,1}
    fwhm::Array{Float64,1}
end

"""
    LaserModel(λlaser::Array{Float64,1},amplitude::Array{Float64,1},fwhm::Array{Float64,1})

Constructor of the model of the laser illumination.
It consist on:
* `λlaser` an array of the wavelengths of the lasers
* `amplitude` an array of the amplitude of the maximum of each Gaussian spots
* `fwhm` an array of the full width at half maximum of the Gaussians
"""
function LaserModel(λlaser::Array{Float64,1},amplitude::Array{Float64,1},fwhm::Array{Float64,1})
    nλ = length(λlaser);
    @assert length(amplitude)==nλ "amplitude vector does not have the right size"
    @assert length(fwhm)==nλ "fwhm vector does not have the right size"
    LaserModel(nλ ,λlaser,amplitude,fwhm);
end

"""
    UpdateLaserModel(self::LaserModel,A::Array{Float64,1},fwhm::Array{Float64,1})

Update the parameters of the laser model
* `self` :  LaserModel object
* `A` : 1D  array of amplitudes of the maximum of each Gaussian spot
* `fwhm` 1D  array of the full width at half maximum of each Gaussian spot

`A` and `fwhm` must have lenth of `self.nλ`
"""
function UpdateLaserModel(self::LaserModel,A::Array{Float64,1},fwhm::Array{Float64,1})
    @assert length(A)==self.nλ "amplitude vector does not have the right size"
    @assert length(fwhm)==self.nλ "fwhm vector does not have the right size"
    self.amplitude = A;
    self.fwhm = fwhm;
end

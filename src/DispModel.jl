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
    cx::Vector{Float64} # coefficients of the polynomial along the x axis
    cy::Vector{Float64} # coefficients of the polynomial along the y axis

    function DispModel(λ0, order, cx, cy)
        order ≥ 0               || throw(ArgumentError)
        length(cx) == (order+1) || throw(ArgumentError)
        length(cy) == (order+1) || throw(ArgumentError)
        new(λ0, order, cx, cy)
    end
end

function DispModel(λ0::Float64, order::Int)
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



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
    updateDispModel(::DispModel, C::Matrix{Float64})

Update the coefficients  of the DispModel .
* `C` : array containing the polynomial coefficients.
"""
function updateDispModel(dmodel::DispModel, C::Matrix{Float64})
    size(C) == (2, dmodel.order+1) || throw(DimensionMismatch(
        "coefficients array does not have the right size"))
    dmodel.cx = C[1,:]
    dmodel.cy = C[2,:]
    return self
end

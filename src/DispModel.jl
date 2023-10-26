
"""
    DispModel(λ0::Float64, order::Int, cx::Vector{Float64}, cy::Vector{Float64})

The dispersion model giving the position of a wavelength on the detector
* `λ0` is the reference wavelength
* `order` is the order of the polynomials
* `cx` is a vector of coefficients of the polynomial along the x axis
* `cy` is a vector of coefficients of the polynomial along the y axis
"""
mutable struct DispModel
    λ0    ::Float64 # reference wavelength
    order ::Int     # order of the polynomial
    cx    ::Vector{Float64} # coefficients of the polynomial along the x axis
    cy    ::Vector{Float64} # coefficients of the polynomial along the y axis

    function DispModel(λ0, order, cx, cy)
        length(cx) == order + 1 || throw(DimensionMismatch(
            "coefficients length `cx` does not match the order"))
        length(cy) == order + 1 || throw(DimensionMismatch(
            "coefficients length `cy` does not match the order"))
        new(λ0, order, cx, cy)
    end
end

"""
    DispModel(λ0::Float64, order::Int)

Version with automatic setting of `cx` and `cy` with zeroes.
"""
function DispModel(λ0, order)
    cx = zeros(order + 1)
    cy = zeros(order + 1)
    DispModel(λ0, order, cx, cy)
end

"""
    (::DispModel)(λ::Float64) -> NTuple{2,Float64}

compute the position `(x,y)` of the wavelength `λ` according to the dispersion law `DispModel`.
"""
function (self::DispModel)(λ::Float64) ::NTuple{2,Float64}

    λpowers = ( (λ-self.λ0) / self.λ0 ).^(1:self.order)
    x = self.cx[1] + sum(self.cx[2:end] .* λpowers)
    y = self.cy[1] + sum(self.cy[2:end] .* λpowers)
    (x, y)
end

"""
    updateDispModel(::DispModel, cx::Vector{Float64}, cy::Vector{Float64})

Update the coefficients of the DispModel.
"""
function updateDispModel(dmodel::DispModel, cx::Vector{Float64}, cy::Vector{Float64})
    length(cx) == dmodel.order + 1 || throw(DimensionMismatch(
        "coefficients vector `cx` does not have the right length"))
    length(cy) == dmodel.order + 1 || throw(DimensionMismatch(
        "coefficients vector `cy` does not have the right length"))
    dmodel.cx = cx
    dmodel.cy = cy
end

mutable struct ProfileModel
    λ0    ::Float64   # reference wavelength
    order ::Int       # order of the polynomial
    cx    ::Vector{Float64} # coefficients of the polynomial along the x axis
    cλ    ::Vector{Float64} # coefficients of the polynomial along the wavelength axis
    
    function ProfileModel(λ0, order, cx, cλ)
        length(cx) == order + 1 || throw(DimensionMismatch(
            "coefficients length `cx` does not match the order"))
        length(cλ) == order + 1 || throw(DimensionMismatch(
            "coefficients length `cλ` does not match the order"))
        new(λ0, order, cx, cλ)
    end
end

function ProfileModel(λ0::Float64, order::Int)
    cx = fill(NaN, order + 1)
    cλ = fill(NaN, order + 1)
    ProfileModel(λ0, order, cx, cλ)
end

function (self::ProfileModel)(λ::Float64)
    x = self.cx[1];
    w = self.cλ[1];
    @inbounds for o in 1:self.order
        λpo = (( λ - self.λ0)/self.λ0 )^(o)
        x += self.cx[o + 1] * λpo;
        w += self.cλ[o + 1] * λpo;
    end
    return (w, x)
end

function (self::ProfileModel)(λ::Float64, z)

    λpowers = ( (λ-self.λ0) / self.λ0 ).^(1:self.order)
    x = self.cx[1] + sum(self.cx[2:end] .* λpowers)
    w = self.cλ[1] + sum(self.cλ[2:end] .* λpowers)
    
    return (w, (x - z)^2)
end

function updateProfileModel(pmodel::ProfileModel, cx::Vector{Float64}, cλ::Vector{Float64})
    length(cx) == pmodel.order + 1 || throw(DimensionMismatch(
        "coefficients vector `cx` does not have the right length"))
    length(cλ) == pmodel.order + 1 || throw(DimensionMismatch(
        "coefficients vector `cλ` does not have the right length"))
    pmodel.cx = cx
    pmodel.cλ = cλ
end

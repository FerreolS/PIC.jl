mutable struct ProfileModel
    λ0::Float64   # reference wavelength
    order::Int  # order of the polynomial
    cλ::Vector{Float64} # coefficients of the polynomial along the wavelength axis
    cy::Vector{Float64} # coefficients of the polynomial along the y axis
    function ProfileModel(λ0,order,cλ,cy)
        order ≥ 0        || throw(ArgumentError)
        length(cλ) == (order+1) || throw(ArgumentError)
        length(cy) == (order+1) || throw(ArgumentError)
        new(λ0,order,cλ,cy)
    end
end

function ProfileModel(λ0::Float64, order::Int)
    cλ = zeros(order+1)
    cλ[1]=1
    cy = zeros(order+1)
    cy[1]=1
    ProfileModel(λ0,order,cλ,cy)
end

function ProfileModel(λ0::Float64, C::Matrix{Float64})
	order = size(C,2)-1
    cλ = C[1,:]
    cy = C[2,:]
    ProfileModel(λ0,order,cλ,cy)
end

function (self::ProfileModel)(λ::Float64, x) ::NTuple{2,Float64}
    λpo = ((λ-self.λ0)/self.λ0).^(1:self.order)
    w = self.cλ[1] + sum(self.cλ[2:end]  .* λpo)
    y = self.cy[1] + sum(self.cy[2:end] .* λpo)
    return (w, (y - x)^2)
end

function UpdateProfileModel(self::ProfileModel, C::Matrix{Float64}) ::ProfileModel
    size(C) == (2,self.order+1) || error("coefficients size does not match the order")
    self.cλ = C[1,:]
    self.cy = C[2,:]
    return self
end

struct Profile_LKL{T<:Real,A<:AbstractMatrix{T},B<:AbstractMatrix{T}}
    model::ProfileModel
    data::A
    weight::B
    λMap::Matrix{T}
    bbox::BoundingBox{Int}
    amplitude::Vector{T}
    function Profile_LKL{T,A,B}(model,data,weight,λMap,bbox,amplitude) where {T,A,B}
        size(data) == size(weight) || throw(ArgumentError)
        size(data) == size(λMap)   || throw(ArgumentError)
        size(data) == size(bbox)   || throw(ArgumentError)
        size(data,2)+1 == length(amplitude) || throw(ArgumentError)
        new{T,A,B}(model, data, weight, λMap, bbox, amplitude)
    end
end

function Profile_LKL(
    model::ProfileModel, data::A, weight::B, λMap::Matrix{T}, bbox::BoundingBox{Int}
) where {T<:Real,A<:AbstractMatrix{T},B<:AbstractMatrix{T}}
    amplitude = zeros(T, size(data,2)+1)
    Profile_LKL{T,A,B}(model, data, weight, λMap, bbox, amplitude)
end

function (self::Profile_LKL)(coefs::Matrix{T}) ::T where {T<:Real}
    UpdateProfileModel(self.model, coefs)
    p = @. GaussianModel2(self.model(self.λMap,($(axes(self.bbox,1)))))
    profile = p ./ sum(p,dims=1)
    amp = Zygote.@ignore  updateAmplitudeAndBackground(profile,self.data,self.weight)
    Zygote.@ignore self.amplitude .= amp[:]
    return (sum(abs2,@. self.weight * (self.data - amp[1] - $(reshape(amp[2:end],1,:)) * profile)))
end

mutable struct ProfileModel
    λ0::Float64   # reference wavelength
    order::Int  # order of the polynomial
    cλ::Vector{Float64} # coefficients of the polynomial along the wavelength axis
    cy::Vector{Float64} # coefficients of the polynomial along the y axis
    function ProfileModel(λ0,order,cλ,cy)
        order ≥ 0               || throw(ArgumentError)
        length(cλ) == (order+1) || throw(ArgumentError)
        length(cy) == (order+1) || throw(ArgumentError)
        new(λ0, order, cλ, cy)
    end
end

function ProfileModel(λ0::Float64, order::Int)
    cλ = zeros(order+1)
    cy = zeros(order+1)
    cλ[1] = 1
    cy[1] = 1
    ProfileModel(λ0, order, cλ, cy)
end

function ProfileModel(λ0::Float64, coefs::Vector{Float64})
    order = Int(length(coefs) / 2) - 1
    cλ = coefs[1 : (order+1)]
    cy = coefs[(order+2) : end]
    ProfileModel(λ0, order, cλ, cy)
end

function (self::ProfileModel)(λ::Float64, x::Int) ::NTuple{2,Float64}
    λpo = ((λ-self.λ0)/self.λ0).^(1:self.order)
    w = self.cλ[1] + sum(self.cλ[2:end] .* λpo)
    y = self.cy[1] + sum(self.cy[2:end] .* λpo)
    return (w, (y - x)^2)
end

function updateProfileModel!(self::ProfileModel, coefs::Vector{Float64}) ::Nothing
    length(coefs) == 2 * (self.order + 1) || throw(ArgumentError)
    self.cλ = coefs[1 : (self.order+1)]
    self.cy = coefs[(self.order+2) : end]
    nothing
end

struct Profile_LKL{D<:AbstractMatrix{<:Real},W<:AbstractMatrix{<:Real}}
    bbox::BoundingBox{Int}
    profile_model::ProfileModel
    data::D
    weights::W
    λMap::Matrix{Float64}
    amplitude::Vector{Float64}
    function Profile_LKL{D,W}(
        bbox, profile_model, data, weights, λMap, amplitude
    ) where {D,W}
        size(data) == size(weights)         || throw(ArgumentError)
        size(data) == size(λMap)            || throw(ArgumentError)
        size(data) == size(bbox)            || throw(ArgumentError)
        size(data,2)+1 == length(amplitude) || throw(ArgumentError)
        new{D,W}(bbox, profile_model, data, weights, λMap, amplitude)
    end
end

function Profile_LKL(
    bbox::BoundingBox{Int}, profile_model::ProfileModel,
    data::D, weight::W, λMap::Matrix{Float64}
) where {D<:AbstractMatrix{<:Real},W<:AbstractMatrix{<:Real}}
    amplitude = zeros(Float64, size(data,2)+1)
    Profile_LKL{D,W}(bbox, profile_model, data, weight, λMap, amplitude)
end

function (self::Profile_LKL)(coefs::Vector{Float64}) ::Float64
    updateProfileModel!(self.profile_model, coefs)
    p = @. GaussianModel2(self.profile_model(self.λMap,($(axes(self.bbox,1)))))
    profile = p ./ sum(p; dims=1)
    amp = Zygote.@ignore updateAmplitudeAndBackground!(profile, self.data, self.weights)
    Zygote.@ignore self.amplitude .= amp[:]
    return (sum(abs2,@. self.weights * (self.data - amp[1] - $(reshape(amp[2:end],1,:)) * profile)))
end

function updateAmplitudeAndBackground!(
    profile, data::D, weights::W
) ::Vector{Float64} where {D<:AbstractMatrix{<:Real},W<:AbstractMatrix{<:Real}}
    
    c = @. profile *  weights
    b = @. profile * data * weights
    a = @. profile^2 * weights
    a = reshape(sum(a; dims=1), Val(1))
    b = reshape(sum(b; dims=1), Val(1))
    c = reshape(sum(c; dims=1), Val(1))
    za = (a .== 0) .|| (b .<= 0)
    if any(za)
        a[za] .= 1
        b[za] .= 0
        c[za] .= 0
    end

    N = length(a)
    A = Matrix{Float64}(undef,N+1,N+1)
    A[1,1] = sum(weights)
    A[1,2:end] .= A[2:end,1] .= c[:]
    A[2:end,2:end] .= diagm(a)

    b = vcat(sum(data .* weights), b[:])

    return inv(A) * b
end

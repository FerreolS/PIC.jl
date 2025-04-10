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
    cλ[1] = 1
    cy = zeros(order+1)
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

struct Profile_LKL{T<:Real,MD<:AbstractMatrix{T},MW<:AbstractMatrix{T}}
    bbox::BoundingBox{Int}
    profile_model::ProfileModel
    data::MD
    weights::MW
    λMap::Matrix{T}
    amplitude::Vector{T}
    function Profile_LKL{T,MD,MW}(
        bbox, profile_model, data, weights, λMap, amplitude
    ) where {T,MD,MW}
        size(data) == size(weights)         || throw(ArgumentError)
        size(data) == size(λMap)            || throw(ArgumentError)
        size(data) == size(bbox)            || throw(ArgumentError)
        size(data,2)+1 == length(amplitude) || throw(ArgumentError)
        new{T,MD,MW}(bbox, profile_model, data, weights, λMap, amplitude)
    end
end

function Profile_LKL(
    bbox::BoundingBox{Int}, profile_model::ProfileModel, data::MD, weight::MW, λMap::Matrix{T}
) where {T<:Real,MD<:AbstractMatrix{T},MW<:AbstractMatrix{T}}
    amplitude = zeros(T, size(data,2)+1)
    Profile_LKL{T,MD,MW}(bbox, profile_model, data, weight, λMap, amplitude)
end

function (self::Profile_LKL)(coefs::Vector{T}) ::T where {T<:Real}
    updateProfileModel!(self.profile_model, coefs)
    p = @. GaussianModel2(self.profile_model(self.λMap,($(axes(self.bbox,1)))))
    profile = p ./ sum(p; dims=1)
    amp = Zygote.@ignore updateAmplitudeAndBackground!(profile, self.data, self.weights)
    Zygote.@ignore self.amplitude .= amp[:]
    return (sum(abs2,@. self.weights * (self.data - amp[1] - $(reshape(amp[2:end],1,:)) * profile)))
end

function updateAmplitudeAndBackground!(
    profile,data::MA,weights::MB
) ::Vector{T} where {T<:AbstractFloat,MA<:AbstractMatrix{T},MB<:AbstractMatrix{T}}
    
    c = @. profile *  weights
    b = @. profile * data * weights
    a = @. profile^2 * weights
    a = sum(a,dims=1)[:]
    b = sum(b,dims=1)[:]
    c = sum(c,dims=1)[:]
    za = (a .== T(0)).||(b.<=T(0))
    if any(za)
        a[za] .=T(1)
        b[za] .=T(0)
        c[za] .=T(0)
    end
    

    N = length(a)
    A = Matrix{T}(undef,N+1,N+1)
    A[1,1] = sum(weights)
    A[1,2:end] .= A[2:end,1] .= c[:]
    A[2:end,2:end] .= diagm(a)

    b =  vcat(sum(data .* weights),b[:])

    return  inv(A)*b
end

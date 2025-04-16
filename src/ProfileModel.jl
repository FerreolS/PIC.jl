mutable struct ProfileModel
    λref::Float64   # reference wavelength
    order::Int  # order of the polynomial
    cλs::Vector{Float64} # coefficients of the polynomial along the wavelength axis
    cxs::Vector{Float64} # coefficients of the polynomial along the x axis
    function ProfileModel(λref,order,cλs,cxs)
        order ≥ 0               || throw(ArgumentError)
        length(cλs) == (order+1) || throw(ArgumentError)
        length(cxs) == (order+1) || throw(ArgumentError)
        new(λref, order, cλs, cxs)
    end
end

function ProfileModel(λref::Float64, order::Int)
    cλs = zeros(order+1)
    cxs = zeros(order+1)
    cλs[1] = 1
    cxs[1] = 1
    ProfileModel(λref, order, cλs, cxs)
end

function ProfileModel(λref::Float64, coefs::Vector{Float64})
    order = Int(length(coefs) / 2) - 1
    cλs = coefs[1 : (order+1)]
    cxs = coefs[(order+2) : end]
    ProfileModel(λref, order, cλs, cxs)
end

function (self::ProfileModel)(λ::Float64, x::Int) ::NTuple{2,Float64} # [?, pix]
    λpo = ((λ-self.λref)/self.λref).^(1:self.order)
    w = self.cλs[1] + sum(self.cλs[2:end] .* λpo)
    gaussian_cxs = self.cxs[1] + sum(self.cxs[2:end] .* λpo) # [pix coord]
    dist_to_gaussian_cxs = (gaussian_cxs - x)^2             # [pix]
    return (w, dist_to_gaussian_cxs)
end

function updateProfileModel!(
    self::ProfileModel, cλs::Vector{Float64}, cxs::Vector{Float64}
) ::Nothing
    length(cλs) == (self.order+1) || throw(ArgumentError)
    length(cxs) == (self.order+1) || throw(ArgumentError)
    self.cλs = cλs
    self.cxs = cxs
    nothing
end

struct Profile_LKL{D<:AbstractMatrix{<:Real},W<:AbstractMatrix{<:Real}}
    bbox::BoundingBox{Int}
    profile_model::ProfileModel
    data::D
    weights::W
    λMap::AbstractMatrix{Float64}
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
    data::D, weight::W, λMap::AbstractMatrix{Float64}
) where {D<:AbstractMatrix{<:Real},W<:AbstractMatrix{<:Real}}
    amplitude = zeros(Float64, size(data,2)+1)
    Profile_LKL{D,W}(bbox, profile_model, data, weight, λMap, amplitude)
end

function encode_profile_lkl_fitvars(cλs::Vector{Float64}, cxs::Vector{Float64}) ::Vector{Float64}
    fitvars = [ cλs ; cxs ]
end

function decode_profile_lkl_fitvars(fitvars::Vector{Float64}) ::NTuple{2,Vector{Float64}}
    half = Int(length(fitvars)/2)
    cλs = fitvars[1:half]
    cxs = fitvars[half+1:end]
    (cλs,cxs)
end

function (self::Profile_LKL)(fitvars::Vector{Float64}) ::Float64
    (cλs,cxs) = decode_profile_lkl_fitvars(fitvars)
    updateProfileModel!(self.profile_model, cλs, cxs)
    lens_rx = axes(self.bbox, 1)
    p = @. GaussianModel2(self.profile_model(self.λMap, $lens_rx))
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

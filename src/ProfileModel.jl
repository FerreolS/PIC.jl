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

function compute_lamp_fwhm_and_center_x(
    order::Int, λref::Float64, cλs::Vector{Float64}, cxs::Vector{Float64}, λ::Float64, x::Int
) ::NTuple{2,Float64}
    λpo = ((λ - λref) / λref).^(1:order)
    fwhm = cλs[1] + sum(cλs[2:end] .* λpo)
    center_x = cxs[1] + sum(cxs[2:end] .* λpo)
    sq_dist_to_center_x = (center_x - x)^2
    (fwhm, sq_dist_to_center_x)
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

struct Lamp_LKL{D<:AbstractMatrix{<:Real},W<:AbstractMatrix{<:Real}}
    order::Int
    λref::Float64
    bbox::BoundingBox{Int}
    lasers_pixels_λs::AbstractMatrix{Float64}
    data::D
    weights::W
    function Lamp_LKL{D,W}(order, λref, bbox, lasers_pixels_λs, data, weights) where {D,W}
        size(bbox) == size(lasers_pixels_λs) || throw(ArgumentError)
        size(bbox) == size(data)             || throw(ArgumentError)
        size(bbox) == size(weights)          || throw(ArgumentError)
        new{D,W}(order, λref, bbox, lasers_pixels_λs, data, weights)
    end
end

function Lamp_LKL(
    order::Int, λref::Float64, bbox::BoundingBox{Int},
    lasers_pixels_λs::AbstractMatrix{Float64}, data::D, weights::W
) where {D<:AbstractMatrix{<:Real},W<:AbstractMatrix{<:Real}}
    Lamp_LKL{D,W}(order, λref, bbox, lasers_pixels_λs, data, weights)
end

function encode_lamp_lkl_fitvars(cλs::Vector{Float64}, cxs::Vector{Float64}) ::Vector{Float64}
    fitvars = [ cλs ; cxs ]
end

function decode_lamp_lkl_fitvars(fitvars::Vector{Float64}) ::NTuple{2,Vector{Float64}}
    half = length(fitvars) ÷ 2
    cλs = fitvars[1:half]
    cxs = fitvars[half+1:end]
    (cλs, cxs)
end

function (self::Lamp_LKL)(fitvars::Vector{Float64}) ::Float64
    (cλs, cxs) = decode_lamp_lkl_fitvars(fitvars)
#    updateProfileModel!(self.profile_model, cλs, cxs)
    bbox_rx = axes(self.bbox, 1)
    
    lamp_image = [
        GaussianModel2(compute_lamp_fwhm_and_center_x(
            self.order, self.λref, cλs, cxs, self.lasers_pixels_λs[x,y], bbox_rx[x]))
        for x in 1:size(self.bbox,1), y in 1:size(self.bbox,2) ]
    
#    p = @. GaussianModel2(self.profile_model(self.lasers_pixels_λs, bbox_rx))
#    println(size(p))
#    error()
    lamp_image_norm = lamp_image ./ sum(lamp_image; dims=1)
    amp = Zygote.@ignore updateAmplitudeAndBackground!(lamp_image_norm, self.data, self.weights)
#    Zygote.@ignore self.amplitude .= amp[:]
    return (sum(abs2,@. self.weights * (self.data - amp[1] - $(reshape(amp[2:end],1,:)) * lamp_image_norm)))
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

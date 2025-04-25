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
    bbox_rx = axes(self.bbox, 1)
    
    lamp_image = [
        GaussianModel2(compute_lamp_fwhm_and_center_x(
            self.order, self.λref, cλs, cxs, self.lasers_pixels_λs[x,y], bbox_rx[x]))
        for x in 1:size(self.bbox,1), y in 1:size(self.bbox,2) ]
    
    lamp_image_norm = lamp_image ./ sum(lamp_image; dims=1)

    amp = updateAmplitudeAndBackground!(lamp_image_norm, self.data, self.weights)
    return (sum(abs2,@. self.weights * (self.data - amp[1] - $(reshape(amp[2:end],1,:)) * lamp_image_norm)))
end

function updateAmplitudeAndBackground!(
    profile, data::D, weights::W
) ::Vector{Float64} where {D<:AbstractMatrix{<:Real},W<:AbstractMatrix{<:Real}}
    
    c = @. profile * weights
    b = @. profile * data * weights
    a = @. profile^2 * weights
    
    va = sum(a; dims=1)[:]
    vb = sum(b; dims=1)[:]
    vc = sum(c; dims=1)[:]
    
    za = (va .== 0) .|| (vb .<= 0)

    va2 = map(i -> za[i] ? 1 : va[i], eachindex(va))
    vb2 = map(i -> za[i] ? 0 : vb[i], eachindex(vb))
    vc2 = map(i -> za[i] ? 0 : vc[i], eachindex(vc))

    N = length(va2)
    A = Matrix{Float64}(undef,N+1,N+1)
    A = hcat( vcat(sum(weights), vc2), vcat(vc2', diagm(va2)) )

    vb3 = vcat(sum(data .* weights), vb2[:])

    return inv(A) * vb3
end

# Polynomes

function compute_polynome(coefs::Vector{Float64}, x::Float64) ::Float64
    powered_xs = x .^ (0 : length(coefs)-1)
    sum(coefs .* powered_xs)

end

function polynome(coefs::Vector{Float64}) ::Function
    (x::Float64) -> compute_polynome(coefs, x)
end

# Polynomes with reference for x

function compute_polynome_with_reference(ref::Float64, coefs::Vector{Float64},x::Float64) ::Float64
    refed_x = (x - ref) / ref
    compute_polynome(coefs, refed_x)
end

function polynome_with_reference(ref::Float64, coefs::Vector{Float64}) ::Function
    (x::Float64) -> compute_polynome_with_reference(ref, coefs, x)
end

# Gaussians

"""
Compute the value at position sqrt(r) 1D centered Gaussian
* `x`:  squared sampled position
* `fwhm` : full-width at half maximum
"""
function Gaussian(x::T, fwhm::T) where {T<:Real}
    fwhm2sigma = T(1 / (2 * sqrt(2 * log(2))))
    return exp(-x / (2 * (fwhm * fwhm2sigma)^2))
end

# Mean of data & weights

function mean_data_and_weights(
    data::Array{T,N}, weights::Array{T,N}) ::NTuple{2,Matrix{T}} where {T,N}

    2 ≤ N ≤ 3                            || throw(DimensionMismatch())
    size(data) == size(weights)          || throw(DimensionMismatch())
    
    nx = size(data, 1)
    ny = size(data, 2)
    nf = size(data, 3)
    
    nf == 1 && return (reshape(copy(data), nx, ny), reshape(copy(weights), nx, ny))
    
    mean_data    = Matrix{T}(undef, nx, ny)
    mean_weights = similar(mean_data)
    
    for y in 1:ny, x in 1:nx
    
        non_zeros = findall(!iszero, data[x,y,:])
        if isempty(non_zeros)
            mean_data[x,y]    = 0
            mean_weights[x,y] = 0
        else
            mean_data[x,y]    = mean(data[x,y,non_zeros])
            mean_weights[x,y] = length(non_zeros) / sum(1 ./ weights[x,y,non_zeros])
        end
    end
    
    mean_data, mean_weights
end

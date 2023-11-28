# Polynomes

"""
    compute_polynome(coefs::Vector{Float64}, x::Float64) -> Float64

Compute the polynome value for variable `x` and coefficients `coefs`.
"""
function compute_polynome(coefs::Vector{Float64}, x::Float64) ::Float64
    powered_xs = x .^ (0 : length(coefs)-1)
    sum(coefs .* powered_xs)

end

"""
    compute_polynome(coefs::Vector{Float64}) -> Function

Return a function of `x`, computing the polynome with the given `coefs`.

See [`compute_polynome`](@ref).
"""
function polynome(coefs::Vector{Float64}) ::Function
    (x::Float64) -> compute_polynome(coefs, x)
    # Note that the anonymous function seems faster than using Fix1
end

# Polynomes with reference for x

"""
    compute_polynome_with_reference(ref::Float64, coefs::Vector{Float64}, x::Float64) -> Float64

Compute the polynome value for variable `x` and coefficients `coefs` and reference `ref`.

The polynome is computed with variable `(x - ref) / ref` instead of just `x`.

The point is to reduce the size of the variable, when we are fitting a "flat" part of a polynome.
"""
function compute_polynome_with_reference(ref::Float64, coefs::Vector{Float64},x::Float64) ::Float64
    refed_x = (x - ref) / ref
    compute_polynome(coefs, refed_x)
end

"""
    polynome_with_reference(ref::Float64, coefs::Vector{Float64}) -> Function

Return a function of `x`, computing the polynome with given `ref` and `coefs`.

See [`compute_polynome_with_reference`](@ref).
"""
function polynome_with_reference(ref::Float64, coefs::Vector{Float64}) ::Function
    (x::Float64) -> compute_polynome_with_reference(ref, coefs, x)
end

# Gaussians

"""
    Gaussian(x::T, fwhm::T) -> T where {T<:Real}

Compute the value at position sqrt(x2) 1D centered Gaussian
* `x2`:  squared sampled position
* `fwhm` : full-width at half maximum
"""
function Gaussian(x2::T, fwhm::T) where {T<:Real}
    fwhm2sigma = T(1 / (2 * sqrt(2 * log(2))))
    return exp(-x2 / (2 * (fwhm * fwhm2sigma)^2))
end

# Mean of data & weights

"""
    mean_data_and_weights(data::Array{T,N}, weights::Array{T,N}) -> NTuple{2,Matrix{T}}

From a data cube, and a corresponding weights cube, return the "weighted mean" data frame and 
weights frame.

For each pixel, only the data frames with weight > 0 are used. The used weights are then combined.

For cubes of depth 1 or matrices, the function only returns back a copy of the sole data frame and
weights frame.
"""
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


using Zygote
using TwoDimensional

# Dispersion model
struct DispModel
    λ::Float64   # reference wavelength
    order::Int32
end

"""
    (x,y) = (self::DispModel)(cx,cy) 
    compute the position (x,y) given by the dispersion law wavelength DispModel.λ 
    using the polynmial interpolation of order DispModel.order and coefficients (cx,cy) .
"""
function (self::DispModel)(cx::Array{Float64,1},cy::Array{Float64,1})
    x= cx[1];
    y= cy[1];
    for o in 2:self.order
        x += cx[o]  * self.λ^(o);
        y += cy[o]  * self.λ^(o);
    end
    return (x,y)
end

"""
    (x,y) = (self::DispModel)(C) 
    compute the position (x,y) given by the dispersion law wavelength DispModel.λ 
    using the polynmial interpolation of order DispModel.order and coefficients (cx=C[1,:],cy=C[2,:]).
    C is an array of size (2,order)
"""
(self::DispModel)(C::Array{Float64,2}) = self(C[1,:],C[2,:])


function  GaussianModel(A::Float64,fwhm::Float64,x::Float64,y::Float64)
    local fwhm2sigma = 1/ (2*sqrt(2*log(2.)))::Float64
    return A * exp( -(x^2 + y^2) / (2 * (fwhm * fwhm2sigma )^2));
end

function GaussianSpot(data::Array{Float64,2},weight::Array{Float64,2},bbox::BoundingBox{Int},A::Float64,fwhm::Float64,x::Float64,y::Float64)
    for I in CartesianIndices(bbox)
        m[I] = GaussianModel(A, fwhm, I[1] -x, I[2] -y);
    end
    return m
end
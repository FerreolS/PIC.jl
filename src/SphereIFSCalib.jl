
using Zygote
using TwoDimensional

struct DispModel
    λ::Float64   # reference wavelength
    order::Int32
end

function XYpos(self::DispModel,cx::Array{Float64,1},cy::Array{Float64,1})
    x= cx[1];
    y= cy[1];
    for o in 2:self.order
        x += cx[o]  * self.λ^(o);
        y += cy[o]  * self.λ^(o);
    end
    return (x,y)
end

XYpos(self::DispModel,C::Array{Float64,2}) = XYpos(self,C[1,:],C[2,:])

# wavelengths
const λ1 = 987.72e-9# laser 1 
const λ2 = 1123.71e-9# laser 2 
const λ3 = 1309.37e-9# laser 3
const λ4 = 1545.10e-9  # laser 4  
const λ0 = λ1# reference

# polynomial law  order
order = 2

laser1pos = DispModel(λ1-λ0,order)
laser2pos = DispModel(λ2-λ0,order)
laser3pos = DispModel(λ3-λ0,order)
laser4pos = DispModel(λ4-λ0,order)

function  GaussianModel(A::Float64,fwhm::Float64,x::Float64,y::Float64)
    local fwhm2sigma = 1/ (2*sqrt(2*log(2.)))::Float64
    return A * exp( -(x^2 + y^2) / (2 * (fwhm * fwhm2sigma )^2));
end



# Example
f(a,fwhm,x) = GaussianModel(a,fwhm,XYpos(laser1pos,x)...)
Cinit = rand(Float64,2,order);
a0= 1.0;
fwhm0= 1.0;
f(a0,fwhm0,Cinit);
gradient(f,a0,fwhm0,Cinit)

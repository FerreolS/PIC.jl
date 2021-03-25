
# Example
f(a,fwhm,x) = GaussianModel(a,fwhm,XYpos(laser1pos,x)...)
Cinit = rand(Float64,2,order);
a0= 1.0;
fwhm0= 1.0;
f(a0,fwhm0,Cinit);
gradient(f,a0,fwhm0,Cinit)

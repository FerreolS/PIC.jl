
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

# Example
f(a,fwhm,x) = GaussianModel(a,fwhm,laser1pos(x)...)
Cinit = rand(Float64,2,order);
a0= 1.0;
fwhm0= 1.0;
f(a0,fwhm0,Cinit);
gradient(f,a0,fwhm0,Cinit)

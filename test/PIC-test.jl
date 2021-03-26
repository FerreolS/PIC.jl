include("../src/SphereIFSCalib.jl")

using TwoDimensional, Zygote

# wavelengths
const λ1 = 987.72e-9# laser 1 
const λ2 = 1123.71e-9# laser 2 
const λ3 = 1309.37e-9# laser 3
const λ4 = 1545.10e-9  # laser 4  
const λ0 = λ1# reference
const λlaser = [λ1,λ2,λ3,λ4]

# model of the LensletModel
bbox = BoundingBox(xmin=1, ymin=1, xmax=10, ymax=10);
lmod = LensletModel(λ0,λlaser,bbox);

data = rand(100,100)
weight = Float64.(rand(100,100).>0.1);


Cinit = rand(Float64,2,4);
a0= rand(Float64,4);
fwhm0= rand(Float64,4);

likelihood(a,fwhm,C) = GaussianSpotsCost(data,weight,lmod,a,fwhm,C)
cost = likelihood(a0,fwhm0,Cinit)

∇cost = gradient(likelihood,a0,fwhm0,Cinit)
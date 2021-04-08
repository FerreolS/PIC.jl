include("../src/SphereIFSCalib.jl")
#using .SphereIFSCalib: LensletModel,GaussianSpotsCost,GaussianSpotsModel,BuildLensletImage
using TwoDimensional, Zygote,StatsBase

# wavelengths
const λ1 = 987.72e-9# laser 1 
const λ2 = 1123.71e-9# laser 2 
const λ3 = 1309.37e-9# laser 3
const λ4 = 1545.10e-9  # laser 4  
const λlaser = [λ1,λ2,λ3,λ4]
const λ0 = mean(λlaser)# reference

# model of the LensletModel
bbox = BoundingBox(xmin=1, ymin=1, xmax=10, ymax=50);

data = rand(100,100)
weight = Float64.(rand(100,100).>0.1);


Cinit = rand(Float64,2,4);
a0= rand(Float64,4);
fwhm0= rand(Float64,4)*10;
laser =  LaserModel(λlaser,a0,fwhm0)
lmod = LensletModel(λ0,laser,bbox);


likelihood(a,fwhm,C) = GaussianSpotsCost(data,weight,lmod,laser,a,fwhm,C)
cost = likelihood(a0,fwhm0,Cinit)

∇cost = gradient(likelihood,a0,fwhm0,Cinit)

m = rand(100,100);
function  g(a::Array{Float64,1},fwhm::Array{Float64,1},C::Array{Float64,2}) 
    return sum(GaussianSpotsModel(lmod,laser,a,fwhm,C))
end

Cinit = [[ 6.2 0 0 0 ]; [20 3e7 0 0]]
UpdateDispModel(lmod.dmodel,Cinit)
limage = LensletLaserImage(lmod,laser)
#heatmap(limage)
data= limage ;
function  g1(a::Array{Float64,1},fwhm::Array{Float64,1},C::Array{Float64,2}) 
    UpdateDispModel(lmod.dmodel, C);
    UpdateLaserModel(laser,a,fwhm);
    return sum((data .- LensletLaserImage(lmod,laser)).^2)
end
gradient(g1,a0,fwhm0,Cinit)
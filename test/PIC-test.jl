using PIC
using Plots, StatsBase, TwoDimensional,OptimPackNextGen

# wavelengths
λ1 = 987.72e-9# laser 1
λ2 = 1123.71e-9# laser 2
λ3 = 1309.37e-9# laser 3
λ4 = 1545.10e-9  # laser 4
λlaser = [λ1,λ2,λ3,λ4]
λ0 = mean(λlaser)# reference

# model of the LensletModel
bbox = BoundingBox(xmin=1, ymin=1, xmax=10, ymax=50);

data = rand(100,100)
weight = Float64.(rand(100,100).>0.1);


Cinit = rand(Float64,2,4);
a0= rand(Float64,4);
fwhm0= rand(Float64,4)*10;
laser =  LaserModel(λlaser,a0,fwhm0)
lmod = LensletModel(λ0,laser.nλ-1,bbox);



Cinit = [[ 6.2 10 0 0 ]; [20 100 0 0]]
UpdateDispModel(lmod.dmodel,Cinit)
limage = LensletLaserImage(lmod,laser)
heatmap(limage)

ldata = view(data, lmod.bbox);
lweight = view(weight, lmod.bbox);
lkl = LikelihoodIFS(lmod,laser, ldata,lweight);
cost(x::Vector{Float64}) = lkl(x)::Float64

xinit = vcat([a0[:],fwhm0[:],Cinit[:]]...);
xopt = vmlmb(cost, xinit; verb=50, ftol=(0.0,0),gtol = (0.0,0));

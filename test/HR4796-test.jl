# Copyright (c) 2021 Ferréol Soulez
#
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

using PIC

using Plots, StatsBase,Statistics, TwoDimensional,OptimPackNextGen, DelimitedFiles
using FITSIO, ProgressMeter
# wavelengths
λ1 = 987.72e-9# laser 1
λ2 = 1123.71e-9# laser 2
λ3 = 1309.37e-9# laser 3
λ4 = 1545.10e-9  # laser 4
λlaser = [λ1,λ2,λ3]
λ0 = mean(λlaser)# reference
wavelengthrange = LinRange(850e-9,1600e-9,50); # coarse wavelength range of the instrument

coeffx = readdlm("/Users/ferreol/Data/SPHERE/HR_4796-HD_95086/coef_pol_x.txt", header = false)
cx0 = coeffx[:,1] .+ 1025;
mcx1 = median(coeffx[:,2])*λ0*1e6;
mcx2 = median(coeffx[:,3])*(λ0*1e6)^2;

coeffy = readdlm("/Users/ferreol/Data/SPHERE/HR_4796-HD_95086/coef_pol_y.txt", header = false)
cy0 = coeffy[:,1].+ 1025;
mcy1 = median(coeffy[:,2])*λ0*1e6;
mcy2 = median(coeffy[:,3])*(λ0*1e6)^2;

lensletnumber= length(cx0)



lampData =  read(FITS("/Users/ferreol/Data/SPHERE/HR_4796-HD_95086/IFS_calib_spec_corrected.fits")[1])
laserData =  read(FITS("/Users/ferreol/Data/SPHERE/HR_4796-HD_95086/IFS_calib_wave_corrected.fits")[1])
badpix = Float64.(read(FITS("/Users/ferreol/Data/SPHERE/HR_4796-HD_95086/IFS_BP_corrected.fits")[1]))
ainit = [1000. , 600. , 200.];
fwhminit = [2.0, 2.0 , 2.0];
laser =  LaserModel(λlaser,ainit,fwhminit);

largeur = 4;
hauteur = 44;

nb_fit = 5 #test sur quelques lenslets

dxmin = 2;
dxmax = 2;
dymin = 21;
dymax = 18;

valid = ((cx0 .- dxmin).>0) .&  ((cx0 .+ dxmax).<2048) .&  ((cy0 .- dymin).>0) .&  ((cy0 .+ dymax).<2048);

lenslettab = Array{Union{LensletModel,Missing}}(missing,lensletnumber);
atab = Array{Union{Float64,Missing}}(missing,3,lensletnumber);
fwhmtab = Array{Union{Float64,Missing}}(missing,3,lensletnumber);
ctab = Array{Union{Float64,Missing}}(missing,2,3,lensletnumber);
p = Progress(lensletnumber; showspeed=true)
Threads.@threads for i in findall(valid)
    bbox = round(Int, BoundingBox(cx0[i,1]-dxmin, cx0[i,1]+dxmax, cy0[i,1]-dymin, cy0[i,1]+dymax));

    lenslettab[i] = LensletModel(λ0,laser.nλ-1, bbox);
    Cinit= [ [cx0[i,1] mcx1 mcx2]; [cy0[i,1] mcy1 mcy2] ];
    xinit = vcat([ainit[:],fwhminit[:],Cinit[:]]...);
    laserDataView = view(laserData, bbox);
    badpixview = view(badpix,bbox)
    lkl = LikelihoodIFS(lenslettab[i],deepcopy(laser), laserDataView,badpixview);
    cost(x::Vector{Float64}) = lkl(x)
    try
        xopt = vmlmb(cost, xinit; verb=false,ftol = (0.0,1e-4),maxeval=500);
        (aopt,fwhmopt,copt) = (xopt[1:(laser.nλ)],xopt[(laser.nλ+1):(2*laser.nλ)],reshape(xopt[(2*laser.nλ+1):(4*laser.nλ)],2,:));
        atab[:,i] = aopt
        fwhmtab[:,i] = fwhmopt
        ctab[:,:,i] = copt
    catch
        continue
    end
    next!(p)
end
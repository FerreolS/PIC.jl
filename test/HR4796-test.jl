# Copyright (c) 2021 Ferréol Soulez
#
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

using Revise, PIC

using Plots, StatsBase,Statistics, DelimitedFiles
using FITSIO
# wavelengths
λ1 = 987.72e-9;# laser 1
λ2 = 1123.71e-9;# laser 2
λ3 = 1309.37e-9;# laser 3
λ4 = 1545.10e-9;  # laser 4
λlaser = [λ1,λ2,λ3];
nλ = length(λlaser);
λ0 = mean(λlaser);# reference
wavelengthrange = LinRange(850e-9,1600e-9,10000); # coarse wavelength range of the instrument

coeffx = readdlm("/Users/ferreol/Data/SPHERE/IFS/HR_4796-HD_95086/coef_pol_x.txt", header = false)
cx0 = coeffx[:,1] .+ 1025;
mcx1 = median(coeffx[:,2])*λ0*1e6;
mcx2 = median(coeffx[:,3])*(λ0*1e6)^2;

coeffy = readdlm("/Users/ferreol/Data/SPHERE/IFS/HR_4796-HD_95086/coef_pol_y.txt", header = false)
cy0 = coeffy[:,1].+ 1025;
mcy1 = median(coeffy[:,2])*λ0*1e6;
mcy2 = median(coeffy[:,3])*(λ0*1e6)^2;


position = hcat(cx0, cy0)
cxinit = [mcx1;mcx2];
cyinit = [mcy1;mcy2];
lensletnumber= length(cx0)



lampData =  read(FITS("/Users/ferreol/Data/SPHERE/IFS/HR_4796-HD_95086/IFS_calib_spec_corrected.fits")[1])
laserData =  read(FITS("/Users/ferreol/Data/SPHERE/IFS/HR_4796-HD_95086/IFS_calib_wave_corrected.fits")[1])
badpix = Float64.(read(FITS("/Users/ferreol/Data/SPHERE/IFS/HR_4796-HD_95086/IFS_BP_corrected.fits")[1]))

fwhminit = [2.3, 2.4 , 2.7];

#largeur = 4;
#hauteur = 44;
dxmin = 2;
dxmax = 2;
dymin = 21;
dymax = 18;
lensletsize = (dxmin, dxmax,dymin,dymax);

valid = ((cx0 .- dxmin).>0) .&  ((cx0 .+ dxmax).<2048) .&  ((cy0 .- dymin).>0) .&  ((cy0 .+ dymax).<2048);


#(lenslettab, atab, fwhmtab,ctab) = fitSpectralLaw(laserData,badpix,λlaser,lensletsize,position,cxinit,cyinit,fwhminit;validlenslets=valid[1:100]);
(lenslettab,  atab, fwhmtab, distweight, λMap) = fitSpectralLawAndProfile(laserData,badpix,lampData,badpix,λlaser,lensletsize,position,cxinit,cyinit,fwhminit,wavelengthrange;validlenslets=valid[1:1000]);

#(lenslettab, distweight, λMap) = fitSpectralLaw(laserData,badpix,λlaser,lensletsize,position,cxinit,cyinit,fwhminit,wavelengthrange;validlenslets=valid);

#=
lenslettab = Array{Union{LensletModel,Missing}}(missing,lensletnumber);
atab = Array{Union{Float64,Missing}}(missing,3,lensletnumber);
fwhmtab = Array{Union{Float64,Missing}}(missing,3,lensletnumber);
ctab = Array{Union{Float64,Missing}}(missing,2,3,lensletnumber);
p = Progress(lensletnumber; showspeed=true)
Threads.@threads for i in findall(valid)
    bbox = round(Int, BoundingBox(cx0[i,1]-dxmin, cx0[i,1]+dxmax, cy0[i,1]-dymin, cy0[i,1]+dymax));

    lenslettab[i] = LensletModel(λ0,nλ-1, bbox);
    Cinit= [ [cx0[i,1] mcx1 mcx2]; [cy0[i,1] mcy1 mcy2] ];
    xinit = vcat([fwhminit[:],Cinit[:]]...);
    laserDataView = view(laserData, bbox);
    badpixview = view(badpix,bbox)
    lkl = LikelihoodIFS(lenslettab[i],λlaser, laserDataView,badpixview);
    cost(x::Vector{Float64}) = lkl(x)
    try
        xopt = vmlmb(cost, xinit; verb=false,ftol = (0.0,1e-8),maxeval=500);
        (fwhmopt,copt) = (xopt[1:(nλ)],reshape(xopt[(nλ+1):(3*nλ)],2,:));
        atab[:,i] = lkl.amplitude;
        fwhmtab[:,i] = fwhmopt
        ctab[:,:,i] = copt
    catch
        continue
    end
    next!(p)
end
ProgressMeter.finish!(p)=#
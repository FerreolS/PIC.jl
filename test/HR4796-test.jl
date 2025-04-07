# Copyright (c) 2021 Ferréol Soulez
#
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

using Revise
using PIC
using StatsBase
using DelimitedFiles
using EasyFITS

ENV["JULIA_DEBUG"] = Main.PIC

# wavelengths
λ1 = 987.72e-9;# laser 1
λ2 = 1123.71e-9;# laser 2
λ3 = 1309.37e-9;# laser 3
λ4 = 1545.10e-9;  # laser 4
λlaser = [λ1,λ2,λ3];
λ0 = mean(λlaser);# reference
wavelengthrange = LinRange(850e-9,1600e-9,10000); # coarse wavelength range of the instrument

coeffx = readdlm("test/HR_4796-HD_95086/coef_pol_x.txt", header = false)
cx0 = coeffx[:,1] .+ 1025;
mcx1 = median(coeffx[:,2])*λ0*1e6;
mcx2 = median(coeffx[:,3])*(λ0*1e6)^2;

coeffy = readdlm("test/HR_4796-HD_95086/coef_pol_y.txt", header = false)
cy0 = coeffy[:,1].+ 1025;
mcy1 = median(coeffy[:,2])*λ0*1e6;
mcy2 = median(coeffy[:,3])*(λ0*1e6)^2;

position = hcat(cx0, cy0)
cxinit = [mcx1;mcx2];
cyinit = [mcy1;mcy2];
lensletnumber= length(cx0)

lampData =  readfits("test/HR_4796-HD_95086/IFS_calib_spec_corrected.fits")
laserData =  readfits(Array{eltype(lampData)}, "test/HR_4796-HD_95086/IFS_calib_wave_corrected.fits")
badpix = readfits(Array{eltype(lampData)}, "test/HR_4796-HD_95086/IFS_BP_corrected.fits")

fwhminit = [2.3, 2.4 , 2.7];

dxmin = 2; # [distance in pixels, from the reference pixel of the lenslet box]
dxmax = 2;
dymin = 21;
dymax = 18;
lensletsize = (dxmin, dxmax,dymin,dymax);

valid = ((cx0 .- dxmin).>0) .&  ((cx0 .+ dxmax).<2048) .&  ((cy0 .- dymin).>0) .&  ((cy0 .+ dymax).<2048);

(lenslettab, laserAmplitude, lampAmplitude, laserfwhm,laserdist, λMap)  = fitSpectralLawAndProfile(laserData,badpix,lampData,badpix,λlaser,λ0,lensletsize,position,cxinit,
    cyinit,fwhminit,wavelengthrange;validlenslets=valid, smalltest=true);

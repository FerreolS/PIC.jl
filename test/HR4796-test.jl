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
lasers_λs = [λ1,λ2,λ3];
λref = mean(lasers_λs);# reference
λrange = LinRange(850e-9,1600e-9,10000); # coarse wavelength range of the instrument

coeffx = readdlm("test/HR_4796-HD_95086/coef_pol_x.txt", header = false)
cx0 = coeffx[:,1] .+ 1025;
mcx1 = median(coeffx[:,2])*λref*1e6;
mcx2 = median(coeffx[:,3])*(λref*1e6)^2;

coeffy = readdlm("test/HR_4796-HD_95086/coef_pol_y.txt", header = false)
cy0 = coeffy[:,1].+ 1025;
mcy1 = median(coeffy[:,2])*λref*1e6;
mcy2 = median(coeffy[:,3])*(λref*1e6)^2;

lenslets_coords = hcat(cx0, cy0)
cxinit = [mcx1;mcx2];
cyinit = [mcy1;mcy2];

lasers_data = readfits("test/HR_4796-HD_95086/IFS_calib_wave_corrected.fits")
lamp_data = readfits(Array{eltype(lasers_data)}, "test/HR_4796-HD_95086/IFS_calib_spec_corrected.fits")
# 1 == good pixel, 0 == bad pixel
good_pixels = readfits(Array{eltype(lasers_data)}, "test/HR_4796-HD_95086/IFS_BP_corrected.fits")

lasers_weights = good_pixels
lamp_weights = good_pixels



fwhminit = [2.3, 2.4 , 2.7];

dxmin = 2; # [distance in pixels, from the reference pixel of the lenslet box]
dxmax = 2;
dymin = 21;
dymax = 18;
lenslet_size = (dxmin, dxmax,dymin,dymax);

valid_lenslets = ((cx0 .- dxmin).>0) .&  ((cx0 .+ dxmax).<2048) .&  ((cy0 .- dymin).>0) .&  ((cy0 .+ dymax).<2048);

(lenslets_models, lasers_amplitudes, lamp_amplitudes, lasers_fwhms, lasers_dists, λmap) = 
    fitSpectralLawAndProfile(
        lasers_data, lasers_weights, lamp_data, lamp_weights, lasers_λs, λref, lenslet_size,
        lenslets_coords, cxinit, cyinit, fwhminit, λrange
        ; valid_lenslets, smalltest=true);


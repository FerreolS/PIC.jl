# Copyright (c) 2021 Ferréol Soulez
#
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

using Revise, PIC
ENV["JULIA_DEBUG"] = Main.PIC

using PyPlot,Statistics, DelimitedFiles
using EasyFITS

using Base: Fix1, Fix2

# wavelengths
λ1 = 987.72e-9;# laser 1
λ2 = 1123.71e-9;# laser 2
λ3 = 1309.37e-9;# laser 3
λ4 = 1545.10e-9;  # laser 4
λlaser = [λ1,λ2,λ3];
nλ = length(λlaser);
λ0 = mean(λlaser);# reference
wavelengthrange = LinRange(850e-9,1600e-9,10000); # coarse wavelength range of the instrument

coeffx = readdlm("coef_pol_x.txt", header = false)
cx0 = coeffx[:,1] .+ 1025;
mcx1 = median(coeffx[:,2])*λ0*1e6;
mcx2 = median(coeffx[:,3])*(λ0*1e6)^2;

coeffy = readdlm("coef_pol_y.txt", header = false)
cy0 = coeffy[:,1].+ 1025;
mcy1 = median(coeffy[:,2])*λ0*1e6;
mcy2 = median(coeffy[:,3])*(λ0*1e6)^2;


position = hcat(cx0, cy0)
cxinit = [mcx1;mcx2];
cyinit = [mcy1;mcy2];
lensletnumber= length(cx0)

laserData = readfits(Array{Float32}, "2580479/reduced_wavespecpos/reduced_SPHER.2020-01-09T12:19:51.620_IFS_WAVE,LAMP_1.650726s_10f_OBS_YJ_IFU.fits")
laserData = mean(laserData; dims=3) |> Fix2(reshape, (2048, 2048))

lampData = readfits(Array{Float32}, "2580479/reduced_wavespecpos/reduced_SPHER.2020-01-09T12:20:56.576_IFS_SPECPOS,LAMP_1.650726s_10f_OBS_YJ_IFU.fits")
lampData = mean(lampData; dims=3) |> Fix2(reshape, (2048, 2048))


lampData = readfits(Array{Float32}, "IFS_calib_spec_corrected.fits")
laserData = readfits(Array{Float32}, "IFS_calib_wave_corrected.fits")
badpix = readfits(Array{Float32}, "IFS_BP_corrected.fits")


laserWeights = readfits(Array{Float32}, "2580479/reduced_wavespecpos/reduced_SPHER.2020-01-09T12:19:51.620_IFS_WAVE,LAMP_1.650726s_10f_OBS_YJ_IFU.fits", ext="weights")
lampWeights = readfits(Array{Float32}, "2580479/reduced_wavespecpos/reduced_SPHER.2020-01-09T12:20:56.576_IFS_SPECPOS,LAMP_1.650726s_10f_OBS_YJ_IFU.fits"; ext="weights")
meanLampData    = zeros(Float32, 2048, 2048)
meanLampWeights = zeros(Float32, 2048, 2048)
meanLaserData    = zeros(Float32, 2048, 2048)
meanLaserWeights = zeros(Float32, 2048, 2048)
for y in 1:2048, x in 1:2048
    goods = findall(!iszero, lampWeights[x,y,:])
    meanLampData[x,y]    = isempty(goods) ? 0 : mean(lampData[x,y,goods])
    meanLampWeights[x,y] = isempty(goods) ? 0 : length(goods) / sum(1 ./ lampWeights[x,y,goods])

    goods = findall(!iszero, laserWeights[x,y,:])
    meanLaserData[x,y]    = isempty(goods) ? 0 : mean(laserData[x,y,goods])
    meanLaserWeights[x,y] = isempty(goods) ? 0 : length(goods) / sum(1 ./ laserWeights[x,y,goods])
end

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
(lenslettab, laserAmplitude, lampAmplitude, laserfwhm,laserdist, λMap)  = fitSpectralLawAndProfile(laserData,laserWeights,lampData,lampWeights,λlaser,lensletsize,position,cxinit,cyinit,fwhminit,wavelengthrange;validlenslets=valid);

(msdLenslettab, msdLaserAmplitude, msdLampAmplitude, msdLaserfwhm, msdLaserdist, msdλMap) =
    fitSpectralLawAndProfile(
        meanLaserData, meanLaserWeights, meanLampData, meanLampWeights, λlaser,
        lensletsize, position, cxinit, cyinit, fwhminit, wavelengthrange
        ; validlenslets=valid);



#function test_results((lenslettab, laserAmplitude, lampAmplitude, laserfwhm,laserdist, λMap))
#    @testset "test_results" begin
#
#    # Amp Laser 1 quantiles(7)
#    goal_laser1_amp_q = [ 1.0e-323, 0.7925838470458985, 0.9788287878036499, 1.0460357666015625,
#                         1.0984456539154053, 1.1487341403961182, 1.2742278671264637 ]
#    test_laser1_amp_q = quantile(laserAmplitude[1,:] |> Fix1(filter, !isnan),
#                                [0.01, 0.10, 0.25, 0.50, 0.75, 0.90, 0.99])
#    @test all(test_laser1_amp_q[i] ≈ goal_laser1_amp_q[i] for i in 1:7)
#
#    # Amp Laser 2 quantiles(7)
#    goal_laser2_amp_q = [ 1.0e-323, 0.6232411861419678, 0.7519862353801727, 0.8382254838943481,
#                         0.9226217567920685, 1.0039421319961548, 1.186977058649063 ]
#    test_laser2_amp_q = quantile(laserAmplitude[2,:] |> Fix1(filter, !isnan),
#                                [0.01, 0.10, 0.25, 0.50, 0.75, 0.90, 0.99])
#    @test all(test_laser2_amp_q[i] ≈ goal_laser2_amp_q[i] for i in 1:7)
#
#    # Amp Laser 3 quantiles(7)
#    goal_laser3_amp_q = [ 1.0e-323, 0.13888812810182571, 0.18341900035738945, 0.22312504798173904,
#                          0.2641787901520729, 0.3051838397979734, 0.3985608237981791 ]
#    test_laser3_amp_q = quantile(laserAmplitude[3,:] |> Fix1(filter, !isnan),
#                                [0.01, 0.10, 0.25, 0.50, 0.75, 0.90, 0.99])
#    @test all(test_laser3_amp_q[i] ≈ goal_laser3_amp_q[i] for i in 1:7)
#
#    # Amp Lamp medians (for each of the 41 values)
#    goal_lamp_amp_meds = [0.8344050777819514, -0.3006410904626612, 0.9473053311509265,
#                          2.973670613245254, 4.936688372440381, 6.163724345677633,
#                          7.004124174001541, 7.998535383086773, 8.944212577419998,
#                          9.370400240701608, 9.499360310282281, 9.80868313908925,
#                          10.410436095016221, 11.002490923963592, 11.497349509082337,
#                          12.045985675532815, 12.63320092812729, 13.15101690734027,
#                          13.342216646654698, 13.148293740461934, 12.960010988822356,
#                          13.02047728929342, 13.384199025611643, 13.968493398861083,
#                          14.467698691385753, 14.681729592824974, 14.82950301225792,
#                          14.987015729441254, 15.124815180789492, 15.164839765807875,
#                          15.154688307677937, 15.164679339310421, 15.165244769956743,
#                          15.14092501068172, 15.237119304816334, 15.286894665768639,
#                          14.882295494538635, 14.119253468339329, 12.379838961380495,
#                          8.050791790215975, 2.6444441873852327]
#    test_lamp_amp_meds = [ median(lampAmplitude[i,:] |> Fix1(filter, !isnan)) for i in 1:41 ]
#    @test all(test_lamp_amp_meds[i] ≈ goal_lamp_amp_meds[i] for i in 1:41)
#
#    # FWHM Laser 1 quantiles(7)
#    goal_laser1_fwhm_q = [0.011445818335227611, 2.1898451916852375, 2.310675540415278,
#                         2.4041311832869208, 2.490676288203164, 2.5900385428749515,
#                         3.3732778166945288]
#    test_laser1_fwhm_q = quantile(laserfwhm[1,:] |> Fix1(filter, !isnan),
#                                [0.01, 0.10, 0.25, 0.50, 0.75, 0.90, 0.99])
#    @test all(test_laser1_fwhm_q[i] ≈ goal_laser1_fwhm_q[i] for i in 1:7)
#
#    # FWHM Laser 2 quantiles(7)
#    goal_laser2_fwhm_q = [6.686993669702557e-5, 2.1584203822840764, 2.2746202880327937,
#                          2.362933269288142, 2.4216033953737597, 2.492183605058115,
#                          3.3420726238450245]
#    test_laser2_fwhm_q = quantile(laserfwhm[2,:] |> Fix1(filter, !isnan),
#                                [0.01, 0.10, 0.25, 0.50, 0.75, 0.90, 0.99])
#    @test all(test_laser2_fwhm_q[i] ≈ goal_laser2_fwhm_q[i] for i in 1:7)
#
#    # FWHM Laser 3 quantiles(7)
#    goal_laser3_fwhm_q = [0.17758348988106548, 2.577537483320824, 2.709659078063296,
#                          2.8598730740852742, 3.089121434680108, 3.580839610565196,
#                          78.22767457103265]
#    test_laser3_fwhm_q = quantile(laserfwhm[3,:] |> Fix1(filter, !isnan),
#                                [0.01, 0.10, 0.25, 0.50, 0.75, 0.90, 0.99])
#    @test all(test_laser3_fwhm_q[i] ≈ goal_laser3_fwhm_q[i] for i in 1:7)
#
#    # λMap quantiles(49)
#    goal_λMap_q = [0.0, 0.0, 0.0, 0.0, 0.0, 9.23957395739574e-7, 9.322832283228323e-7,
#                   9.403840384038403e-7, 9.484848484848485e-7, 9.567356735673566e-7,
#                   9.64986498649865e-7, 9.733123312331234e-7, 9.817131713171318e-7,
#                   9.901890189018902e-7, 9.987398739873988e-7, 1.0073657365736572e-6,
#                   1.016066606660666e-6, 1.0248424842484247e-6, 1.0336933693369335e-6,
#                   1.0426192619261927e-6, 1.0516201620162015e-6, 1.0607710771077107e-6,
#                   1.0699219921992199e-6, 1.079297929792979e-6, 1.0886738673867385e-6,
#                   1.0981998199819981e-6, 1.107725772577258e-6, 1.1174767476747673e-6,
#                   1.1272277227722772e-6, 1.137203720372037e-6, 1.1472547254725472e-6,
#                   1.1573807380738074e-6, 1.1676567656765676e-6, 1.178082808280828e-6,
#                   1.1886588658865884e-6, 1.1993099309930992e-6, 1.21018601860186e-6,
#                   1.2211371137113712e-6, 1.2323132313231321e-6, 1.2435643564356436e-6,
#                   1.2551155115511551e-6, 1.2668166816681668e-6, 1.2785928592859288e-6,
#                   1.2907440744074406e-6, 1.303045304530453e-6, 1.3155715571557154e-6,
#                   1.3283978397839783e-6, 1.3413741374137414e-6, 1.3546504650465046e-6]
#    test_λMap_q = quantile(λMap[:], 0.02:0.02:0.98)
#    @test all(test_λMap_q[i] ≈ goal_λMap_q[i] for i in 1:49)
#
#    end
#end

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
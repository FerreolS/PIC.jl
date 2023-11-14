using Revise
using PIC

using EasyFITS, DelimitedFiles, StatsBase
using Base: Fix1, Fix2

# set your own
cd("/home/antoine/Documents/pictest")

ENV["JULIA_DEBUG"] = Main.PIC

T = Float32

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

fwhminit = [2.3, 2.4 , 2.7];

dxmin = 2;
dxmax = 2;
dymin = 21;
dymax = 18;
lensletsize = (dxmin, dxmax,dymin,dymax);

validlensmap = ((cx0 .- dxmin).>0) .&  ((cx0 .+ dxmax).<2048) .&  ((cy0 .- dymin).>0) .&  ((cy0 .+ dymax).<2048);

laserData = readfits(Array{Float64}, "2580479/reduced_wavespecpos/reduced_SPHER.2020-01-09T12:19:51.620_IFS_WAVE,LAMP_1.650726s_10f_OBS_YJ_IFU.fits")

lampData = readfits(Array{Float64}, "2580479/reduced_wavespecpos/reduced_SPHER.2020-01-09T12:20:56.576_IFS_SPECPOS,LAMP_1.650726s_10f_OBS_YJ_IFU.fits")

#laserData = reshape(mean(laserData; dims=3), (2048, 2048))
#lampData  = reshape(mean(lampData;  dims=3), (2048, 2048))


laserWeights = readfits(Array{Float64}, "2580479/reduced_wavespecpos/reduced_SPHER.2020-01-09T12:19:51.620_IFS_WAVE,LAMP_1.650726s_10f_OBS_YJ_IFU.fits", ext="weights")

lampWeights = readfits(Array{Float64}, "2580479/reduced_wavespecpos/reduced_SPHER.2020-01-09T12:20:56.576_IFS_SPECPOS,LAMP_1.650726s_10f_OBS_YJ_IFU.fits"; ext="weights")



meanLampData    = zeros(Float64, 2048, 2048)
meanLampWeights = zeros(Float64, 2048, 2048)
meanLaserData    = zeros(Float64, 2048, 2048)
meanLaserWeights = zeros(Float64, 2048, 2048)
for y in 1:2048, x in 1:2048
    goods = findall(!iszero, lampWeights[x,y,:])
    meanLampData[x,y]    = isempty(goods) ? 0 : mean(lampData[x,y,goods])
    meanLampWeights[x,y] = isempty(goods) ? 0 : length(goods) / sum(1 ./ lampWeights[x,y,goods])

    goods = findall(!iszero, laserWeights[x,y,:])
    meanLaserData[x,y]    = isempty(goods) ? 0 : mean(laserData[x,y,goods])
    meanLaserWeights[x,y] = isempty(goods) ? 0 : length(goods) / sum(1 ./ laserWeights[x,y,goods])
end



#lampData = readfits("IFS_calib_spec_corrected.fits")
#laserData = readfits("IFS_calib_wave_corrected.fits")
#badpix = readfits(Array{Float64}, "IFS_BP_corrected.fits")


#output = fitSpectralLawAndProfile(
#    meanLaserData, meanLaserWeights, meanLampData, meanLampWeights,
#    λlaser, lensletsize, position, cxinit, cyinit, fwhminit, wavelengthrange
#    ; validlensmap)

#=
msdres = (msdLenslettab, msdLaserAmplitude, msdLampAmplitude, msdLaserfwhm, msdLaserdist, msdλMap) =
    fitSpectralLawAndProfile(
        meanLaserData, meanLaserWeights, meanLampData, meanLampWeights, λlaser,
        lensletsize, position, cxinit, cyinit, fwhminit, wavelengthrange
        ; validlenslets=valid)
=#

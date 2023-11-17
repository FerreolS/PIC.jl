using Revise
using PIC

using EasyFITS, StatsBase

# folder and files

# set your own
cd("/home/antoine/Documents/pictest/4las")

vpm_path = "IFS_BP_corrected.fits"
wavelamps_filepath = "IFS_calib_wave_corrected.fits"
specpos_filepath = "IFS_calib_spec_corrected.fits"

#wavelamps_filepath = "2580479/reduced_wavespecpos/reduced_SPHER.2020-01-09T12:19:51.620_IFS_WAVE,LAMP_1.650726s_10f_OBS_YJ_IFU.fits"
#
#specpos_filepath = "2580479/reduced_wavespecpos/reduced_SPHER.2020-01-09T12:20:56.576_IFS_SPECPOS,LAMP_1.650726s_10f_OBS_YJ_IFU.fits"
#
ENV["JULIA_DEBUG"] = Main.PIC

# wavelengths

wavelamps_λlasers = WAVELAMPS_λLASERS_YJH

λ0 = mean(wavelamps_λlasers)

(cx0, mcxs, cy0, mcys) = get_anthony_cxy(λ0, :YJH)

lenses_positions = vcat(cx0', cy0')
cxinit = mcxs
cyinit = mcys

fwhminit = [2.3, 2.4 , 2.7, 2.8];

valid_lenses_map = compute_first_valid_lenses_map(cx0, cy0, DXMIN, DXMAX, DYMIN, DYMAX)

wavelamps_data = readfits(Array{Float64}, wavelamps_filepath)
wavelamps_weights = readfits(Array{Float64}, vpm_path)

specpos_data = readfits(Array{Float64}, specpos_filepath)
specpos_weights = wavelamps_weights

#wavelamps_data_cube    = readfits(Array{Float64}, wavelamps_filepath)
#wavelamps_weights_cube = readfits(Array{Float64}, wavelamps_filepath; ext="weights")
#
#specpos_data_cube    = readfits(Array{Float64}, specpos_filepath)
#specpos_weights_cube = readfits(Array{Float64}, specpos_filepath; ext="weights")
#
#wavelamps_data    = zeros(Float64, 2048, 2048)
#wavelamps_weights = zeros(Float64, 2048, 2048)
#specpos_data     = zeros(Float64, 2048, 2048)
#specpos_weights  = zeros(Float64, 2048, 2048)
#for y in 1:2048, x in 1:2048
#    valids = findall(!iszero, wavelamps_weights_cube[x,y,:])
#    wavelamps_data[x,y]    = isempty(valids) ? 0 : mean(wavelamps_data_cube[x,y,valids])
#    wavelamps_weights[x,y] = isempty(valids) ? 0 : length(valids) / sum(1 ./ wavelamps_weights_cube[x,y,valids])
#
#    valids = findall(!iszero, specpos_weights_cube[x,y,:])
#    specpos_data[x,y]    = isempty(valids) ? 0 : mean(specpos_data_cube[x,y,valids])
#    specpos_weights[x,y] = isempty(valids) ? 0 : length(valids) / sum(1 ./ specpos_weights_cube[x,y,valids])
#end

output = fit_wavelamps_specpos(
    lenses_positions,
    wavelamps_λlasers, wavelamps_data, wavelamps_weights, fwhminit, cxinit, cyinit,
    specpos_data, specpos_weights
    ; λ0, box_frame = BOX_FRAME, valid_lenses_map)

using Revise
using PIC

using EasyFITS, StatsBase

# folder and files

# set your own
cd("/home/antoine/Documents/pictest/2035882")

#vpm_path = "IFS_BP_corrected.fits"
#wavelamps_filepath = "IFS_calib_wave_corrected.fits"
#specpos_filepath = "IFS_calib_spec_corrected.fits"

wavelamps_filepath = "reduced_wavespecpos/reduced_SPHER.2018-05-07T16:24:18.908_IFS_WAVE,LAMP_1.650726s_10f_OBS_H_IFU.fits"

specpos_filepath = "reduced_wavespecpos/reduced_SPHER.2018-05-07T16:23:02.392_IFS_SPECPOS,LAMP_1.650726s_10f_OBS_H_IFU.fits"
#
ENV["JULIA_DEBUG"] = Main.PIC

# wavelengths

wavelamps_λlasers = WAVELAMPS_λLASERS_YJH

λ0 = mean(wavelamps_λlasers)

(cx0, mcxs, cy0, mcys) = get_anthony_cxy(λ0, :YJH)

lenses_positions = vcat(cx0', cy0')
cxinit = mcxs
cyinit = mcys

valid_lenses_map = compute_first_valid_lenses_map(cx0, cy0, DXMIN, DXMAX, DYMIN, DYMAX)

#wavelamps_data = readfits(Array{Float64}, wavelamps_filepath)
#wavelamps_weights = readfits(Array{Float64}, vpm_path)
#
#specpos_data = readfits(Array{Float64}, specpos_filepath)
#specpos_weights = wavelamps_weights

wavelamps_data_cube    = readfits(Array{Float64}, wavelamps_filepath)
wavelamps_weights_cube = readfits(Array{Float64}, wavelamps_filepath; ext="weights")

specpos_data_cube    = readfits(Array{Float64}, specpos_filepath)
specpos_weights_cube = readfits(Array{Float64}, specpos_filepath; ext="weights")

wavelamps_data, wavelamps_weights = mean_data_and_weights(
    wavelamps_data_cube, wavelamps_weights_cube)

specpos_data, specpos_weights = mean_data_and_weights(
    specpos_data_cube, specpos_weights_cube)
    

output = fit_wavelamps_specpos(
    lenses_positions,
    wavelamps_λlasers, wavelamps_data, wavelamps_weights, WAVELAMPS_INIT_FWHM_YJH, cxinit, cyinit,
    specpos_data, specpos_weights
    ; λ0, box_frame = BOX_FRAME, valid_lenses_map)

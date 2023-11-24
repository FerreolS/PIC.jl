using Revise
using PIC
using EasyFITS
using StatsBase

ENV["JULIA_DEBUG"] = Main.PIC

# set your own
@show cd("/home/antoine/Documents/pictest/2580479/")

# ========= ferreol files ================#

#vpm_path = "IFS_BP_corrected.fits"
#wavelamps_filepath = "IFS_calib_wave_corrected.fits"
#specpos_filepath = "IFS_calib_spec_corrected.fits"
#
#wavelamps_data = readfits(Array{Float64}, wavelamps_filepath)
#wavelamps_weights = readfits(Array{Float64}, vpm_path)
#
#specpos_data = readfits(Array{Float64}, specpos_filepath)
#specpos_weights = wavelamps_weights


# ========= or  ================#
# ========= calib files cleaned by scientific detectors ================#

wavelamps_filepath = "reduced_wavespecpos/reduced_SPHER.2020-01-09T12:19:51.620_IFS_WAVE,LAMP_1.650726s_10f_OBS_YJ_IFU.fits"

specpos_filepath = "reduced_wavespecpos/reduced_SPHER.2020-01-09T12:20:56.576_IFS_SPECPOS,LAMP_1.650726s_10f_OBS_YJ_IFU.fits"

wavelamps_data_cube    = readfits(Array{Float64}, wavelamps_filepath)
wavelamps_weights_cube = readfits(Array{Float64}, wavelamps_filepath; ext="weights")

specpos_data_cube    = readfits(Array{Float64}, specpos_filepath)
specpos_weights_cube = readfits(Array{Float64}, specpos_filepath; ext="weights")

wavelamps_data, wavelamps_weights = mean_data_and_weights(
    wavelamps_data_cube, wavelamps_weights_cube)

specpos_data, specpos_weights = mean_data_and_weights(
    specpos_data_cube, specpos_weights_cube)

# =========  ================#

wavelamps_λlasers = WAVELAMPS_λLASERS_YJ

λ0 = mean(wavelamps_λlasers)

(cx0, mcxs, cy0, mcys) = get_anthony_cxy(λ0, :YJ)

lenses_positions = vcat(cx0', cy0')
cxinit = mcxs
cyinit = mcys

result = fit_wavelamps_specpos(
    lenses_positions,
    wavelamps_λlasers, wavelamps_data, wavelamps_weights, WAVELAMPS_INIT_FWHM_YJ, cxinit, cyinit,
    specpos_data, specpos_weights
    ; λ0)

COSTA Guglielmo 
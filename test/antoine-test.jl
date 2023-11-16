using Revise
using PIC

using EasyFITS, StatsBase

# folder and files

# set your own
cd("/home/antoine/Documents/pictest")

wavelamps_filepath = "2580479/reduced_wavespecpos/reduced_SPHER.2020-01-09T12:19:51.620_IFS_WAVE,LAMP_1.650726s_10f_OBS_YJ_IFU.fits"

specpos_filepath = "2580479/reduced_wavespecpos/reduced_SPHER.2020-01-09T12:20:56.576_IFS_SPECPOS,LAMP_1.650726s_10f_OBS_YJ_IFU.fits"

ENV["JULIA_DEBUG"] = Main.PIC

# wavelengths

wavelamps_λlasers = WAVELAMPS_λLASERS[1:3]

λ0 = mean(wavelamps_λlasers)

(cx0, mcx1, mcx2, cy0, mcy1, mcy2) = get_anthony_cxy(λ0)

lenses_positions = vcat(cx0', cy0')
cxinit = [ mcx1, mcx2]
cyinit = [ mcy1, mcy2]

fwhminit = [2.3, 2.4 , 2.7];

valid_lenses_map = compute_first_valid_lenses_map(cx0, cy0, DXMIN, DXMAX, DYMIN, DYMAX)

wavelamps_data_cube    = readfits(Array{Float64}, wavelamps_filepath)
wavelamps_weights_cube = readfits(Array{Float64}, wavelamps_filepath; ext="weights")

specpos_data_cube    = readfits(Array{Float64}, specpos_filepath)
specpos_weights_cube = readfits(Array{Float64}, specpos_filepath; ext="weights")

wavelamps_data    = zeros(Float64, 2048, 2048)
wavelamps_weights = zeros(Float64, 2048, 2048)
specpos_data     = zeros(Float64, 2048, 2048)
specpos_weights  = zeros(Float64, 2048, 2048)
for y in 1:2048, x in 1:2048
    valids = findall(!iszero, wavelamps_weights_cube[x,y,:])
    wavelamps_data[x,y]    = isempty(valids) ? 0 : mean(wavelamps_data_cube[x,y,valids])
    wavelamps_weights[x,y] = isempty(valids) ? 0 : length(valids) / sum(1 ./ wavelamps_weights_cube[x,y,valids])

    valids = findall(!iszero, specpos_weights_cube[x,y,:])
    specpos_data[x,y]    = isempty(valids) ? 0 : mean(specpos_data_cube[x,y,valids])
    specpos_weights[x,y] = isempty(valids) ? 0 : length(valids) / sum(1 ./ specpos_weights_cube[x,y,valids])
end

output = fit_wavelamps_specpos(
    lenses_positions,
    wavelamps_λlasers, wavelamps_data, wavelamps_weights, fwhminit, cxinit, cyinit,
    specpos_data, specpos_weights
    ; λ0, box_frame = BOX_FRAME, valid_lenses_map)


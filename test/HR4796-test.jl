# Copyright (c) 2021 Ferréol Soulez
#
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

using Revise
using PIC
using StatsBase
using EasyFITS

ENV["JULIA_DEBUG"] = Main.PIC

# wavelengths
nλ = 3
lasers_λs = PIC.LASERS_λS[1:nλ]
lasers_fwhms_init = PIC.LASERS_FWHMS_INIT[1:nλ]
λrange = PIC.λRANGE

lasers_data = readfits("test/HR_4796-HD_95086/IFS_calib_wave_corrected.fits")
lamp_data   = readfits("test/HR_4796-HD_95086/IFS_calib_spec_corrected.fits")

# 1 == good pixel, 0 == bad pixel
good_pixels = readfits("test/HR_4796-HD_95086/IFS_BP_corrected.fits")

lasers_weights = good_pixels
lamp_weights = good_pixels

lens_dx_lower = PIC.LENS_DX_LOWER
lens_dx_upper = PIC.LENS_DX_UPPER
lens_dy_lower = PIC.LENS_DY_LOWER
lens_dy_upper = PIC.LENS_DY_UPPER

valid_lenslets = trues(NLENS)

# testing on a small subset for dev
test_indices = [50, 194, 243, 273, 377, 416, 500, 512, 514, 591, 639, 646, 658, 742, 777, 789, 985, 1083, 1104, 1135, 1162, 1203, 1440, 1441, 1526, 1618, 1628, 1678, 1713, 1813, 1944, 1961, 2091, 2218, 2239, 2262, 2334, 2335, 2336, 2420, 2532, 2568, 2624, 2782, 2797, 2824, 3118, 3219, 3291, 3364, 3423, 3451, 3505, 3566, 3572, 3625, 3722, 3740, 3788, 3794, 3909, 3941, 4028, 4051, 4109, 4443, 4513, 4555, 4617, 4680, 4743, 4876, 4880, 4889, 4970, 4978, 5070, 5099, 5130, 5136, 5138, 5201, 5234, 5408, 5417, 5573, 5588, 5654, 5667, 5703, 5816, 5960, 6014, 6042, 6044, 6087, 6107, 6139, 6259, 6369, 6463, 6494, 6498, 6505, 6522, 6645, 6720, 6820, 6891, 7010, 7035, 7080, 7119, 7156, 7287, 7293, 7373, 7478, 7503, 7517, 7596, 7598, 7751, 7855, 7951, 7974, 7998, 8056, 8157, 8164, 8181, 8227, 8434, 8455, 8496, 8506, 8511, 8725, 8764, 8768, 8868, 9039, 9055, 9125, 9160, 9252, 9345, 9444, 9490, 9493, 9516, 9536, 9540, 9619, 9635, 9786, 9789, 9854, 9879, 9935, 10006, 10024, 10166, 10193, 10334, 10390, 10402, 10464, 10498, 10538, 10643, 10751, 10805, 10862, 10920, 10958, 10985, 11106, 11184, 11343, 11427, 11528, 11667, 11757, 11767, 11921, 11935, 11979, 12023, 12067, 12344, 12354, 12371, 12461, 12511, 12517, 12555, 12659, 12708, 12770, 12852, 12867, 13004, 13238, 13280, 13283, 13343, 13434, 13711, 13731, 13755, 13810, 13819, 13910, 13945, 13983, 13988, 14011, 14171, 14301, 14344, 14388, 14473, 14528, 14566, 14623, 14725, 14925, 15000, 15063, 15277, 15338, 15348, 15395, 15473, 15595, 15657, 15745, 15845, 15987, 16001, 16148, 16156, 16170, 16248, 16254, 16272, 16318, 16443, 16460, 16494, 16631, 16715, 16786, 16823, 16865, 16893, 16925, 16972, 16973, 16982, 17050, 17053, 17074, 17077, 17114, 17121, 17180, 17186, 17394, 17437, 17469, 17478, 17498, 17596, 17615, 17619, 17696, 17863, 17893, 17911, 18125, 18141, 18161, 18171, 18268, 18306, 18313, 18406, 18483, 18576, 18615, 18648, 18691, 18747, 18776, 18802]
valid_lenslets .= false
valid_lenslets[test_indices] .= true

(lenslets_models, lasers_amplitudes, lamp_amplitudes, lasers_fwhms, lasers_dists, λmap) = 
    fitSpectralLawAndProfile(
        lasers_data, lasers_weights, lamp_data, lamp_weights
        ; nλ, lasers_λs, lasers_fwhms_init, λrange, valid_lenslets,
          lens_dx_lower, lens_dx_upper, lens_dy_lower, lens_dy_upper);


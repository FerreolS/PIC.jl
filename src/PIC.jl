# Copyright (c) 2021 Ferréol Soulez
#
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

module PIC

export BOX_FRAME, DXMIN, DXMAX, DYMIN, DYMAX,
       WAVELAMPS_λLASERS_YJH, WAVELAMPS_λLASERS_YJ,
       WAVELAMPS_INIT_FWHM_YJH, WAVELAMPS_INIT_FWHM_YJ,
       mean_data_and_weights,
       get_anthony_cxy,
       compute_first_valid_lenses_map,
       fit_wavelamps_specpos

using Zygote,
      StaticArrays,
      StatsBase,
      LinearAlgebra,
      TwoDimensional,
      ProgressMeter,
      OptimPackNextGen,
      DelimitedFiles

"Default box size that will surround the center of each lens"
const INIT_BOX = BoundingBox{Int}(xmin=-2, xmax=2, ymin=-21, ymax=18)

"Lasers wavelenghts in IFS in Wave,Lamp files for YJH mode"
const WAVELAMPS_λLASERS_YJH = Float64[ 987.72e-9, 1123.71e-9, 1309.37e-9, 1545.10e-9 ]

"Lasers wavelenghts in IFS in Wave,Lamp files for YJ mode"
const WAVELAMPS_λLASERS_YJ  = WAVELAMPS_λLASERS_YJH[1:3]

"Default initial values for the FWHN of lasers in IFS in Wave,Lamp files for YJH mode"
const WAVELAMPS_INIT_FWHM_YJH = Float64[2.6, 2.5 , 2.9, 2.4]

"Default initial values for the FWHN of lasers in IFS in Wave,Lamp files for YJ mode"
const WAVELAMPS_INIT_FWHM_YJ  = Float64[2.38171,  2.30755,  2.39567]

"Wavelength range of the IFS instrument"
const IFS_λRANGE = LinRange(850e-9, 1600e-9, 10_000)

include("utils.jl")
include("init_from_anthony.jl")
include("WaveLampLikelihood.jl")
include("SpecPosLikelihood.jl")
include("FitResult.jl")
include("fit.jl")

end

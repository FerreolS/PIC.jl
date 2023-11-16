# Copyright (c) 2021 Ferréol Soulez
#
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

module PIC

export WAVELAMPS_λLASERS,
       DXMIN, DXMAX, DYMIN, DYMAX,
       BOX_FRAME,
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

const (DXMIN, DXMAX, DYMIN, DYMAX) = (2, 2, 21, 18)

const BOX_FRAME = (DXMIN, DXMAX, DYMIN, DYMAX)

const WAVELAMPS_λLASERS = [ 987.72e-9, 1123.71e-9, 1309.37e-9, 1545.10e-9 ]

const IFS_λRANGE = LinRange(850e-9, 1600e-9, 10_000)

include("utils.jl")
include("init_from_anthony.jl")
include("WaveLampLikelihood.jl")
include("SpecPosLikelihood.jl")
include("fit.jl")

end

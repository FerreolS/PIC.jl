# Copyright (c) 2021 Ferréol Soulez
#
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

module PIC

using Zygote, StaticArrays,StatsBase, LinearAlgebra, EasyFITS, TwoDimensional, ProgressMeter, 
      OptimPackNextGen, Random, DelimitedFiles

export fitSpectralLawAndProfile, exporte, importe, compar

include("LasersModel.jl")
include("ProfileModel.jl")
include("LensletModel.jl")
include("calib.jl")
include("io.jl")

end # module

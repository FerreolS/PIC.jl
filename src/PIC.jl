# Copyright (c) 2021 Ferréol Soulez
#
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

module PIC

using Zygote, StaticArrays, StatsBase, LinearAlgebra, EasyFITS, TwoDimensional, ProgressMeter,
    OptimPackNextGen, Random, DelimitedFiles, ChainRulesCore, WeightedData, ConcreteStructs, DifferentiationInterface,
    Parameters
#Optimisers, ComponentArrays,  DifferentiationInterface

export fitSpectralLawAndProfile, exporte, importe, compar, PICParams

include("lasers.jl")
include("lamps.jl")
include("calib.jl")
include("io.jl")

end # module

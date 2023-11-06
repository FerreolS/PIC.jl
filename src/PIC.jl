# Copyright (c) 2021 Ferréol Soulez
#
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

module PIC

export DispModel,
       LaserModel,
       LensletModel,
       LikelihoodDisp,
       LikelihoodProfile,
       fitSpectralLaw,
       fitSpectralLawAndProfile

using Zygote,
      StaticArrays,
      StatsBase,
      LinearAlgebra,
      TwoDimensional,
      ProgressMeter,
      OptimPackNextGen
      
include("utils.jl")
include("DispModel.jl")
include("ProfileModel.jl")
include("LensletModel.jl")
include("LikelihoodDisp.jl")
include("LikelihoodProfile.jl")
include("WaveLampLikelihood.jl")
include("SpecPosLikelihood.jl")
include("fit.jl")

end

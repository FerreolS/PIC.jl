using Test
using PIC
using InteractiveUtils

using PIC: computeAmplitudeWaveLamp, computeAmplitudeAndBackgroundSpecPos

include("testdata.jl")

@testset "WaveLampLikelihood.jl" begin

    @testset "computeAmplitudeWaveLamp" begin
        (; spots, data, weights, res) = testdata_updateAmplitudeWaveLamp
        test_res = computeAmplitudeWaveLamp(spots, data, weights)
        @test isapprox(test_res, res; atol=0.000_000_001)
    end

end # @testset "WaveLampLikelihood.jl

@testset "SpecPosLikelihood.jl" begin

    @testset "computeAmplitudeAndBackgroundSpecPos" begin
        (; profile, data, weights, res) = testdata_computeAmplitudeAndBackgroundSpecPos
        res_bg, res_amps = res[1], res[2]
        test_res = computeAmplitudeAndBackgroundSpecPos(profile, data, weights)
        test_res_bg, test_res_amps = test_res[1], test_res[2]
        @test isapprox(res_bg, test_res_bg; atol=0.000_000_001)
        @test isapprox(res_amps, test_res_amps; atol=0.000_000_001)
    end

end # @testset "computeAmplitudeAndBackgroundSpecPos"

using FiniteVolumeMethod1D
using Aqua
using JET
using Test

@testset "Quality Assurance" begin
    @testset "Aqua" begin
        Aqua.test_all(FiniteVolumeMethod1D)
    end
    @testset "JET" begin
        JET.test_package(FiniteVolumeMethod1D; target_defined_modules = true)
    end
end

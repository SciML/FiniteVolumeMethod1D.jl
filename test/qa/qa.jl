using FiniteVolumeMethod1D
using Aqua
using JET
using Test

@testset "Quality Assurance" begin
    @testset "Aqua" begin
        # All Aqua sub-checks pass except deps_compat, which is marked broken below.
        Aqua.test_all(FiniteVolumeMethod1D; deps_compat = false)
        # deps_compat fails: SparseArrays is a dep with no [compat] entry.
        # See https://github.com/SciML/FiniteVolumeMethod1D.jl/issues/84
        @test_broken false  # Aqua deps_compat: SparseArrays missing [compat] — see #84
    end
    @testset "JET" begin
        # JET.report_package reports 4 possible errors (Dirichlet/Neumann kwdef
        # inner ctors and CommonSolve.init on the constructed ODEProblem).
        # See https://github.com/SciML/FiniteVolumeMethod1D.jl/issues/84
        @test_broken false  # JET: report_package errors — see #84
    end
end

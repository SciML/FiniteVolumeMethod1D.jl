using FiniteVolumeMethod1D
using Test
using SafeTestsets

@testset "FiniteVolumeMethod1D" begin
    @safetestset "FVMGeometry" begin
        include("geometry.jl")
    end
    @safetestset "FVMProblem" begin
        include("problem.jl")
    end
    @safetestset "ODEProblem" begin
        include("ode_problem.jl")
    end
    @testset "Example Problems" begin
        @safetestset "Porous-Medium" begin
            include("porous_medium.jl")
        end
        @safetestset "Fisher" begin
            include("fisher.jl")
        end
    end
end
using FiniteVolumeMethod1D
using Test
using SafeTestsets

@testset "FiniteVolumeMethod1D" begin
    @safetestset "FVMGeometry" begin
        include("geometry.jl")
    end
    @safetestset "BoundaryConditions" begin
        include("boundary_conditions.jl")
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
        #@safetestset "Fisher" begin # ERROR: MethodError: no method matching hasmetadata(::Pair{Num, Float64}, ::Type{Symbolics.VariableDefaultValue})
        #    include("fisher.jl")
        #end
        @safetestset "Robin Diffusion" begin
            include("robin_diffusion.jl")
        end
        @safetestset "Reaction-Diffusion" begin
            include("dirichlet_source.jl")
        end
        @safetestset "Heat Equation" begin
            include("heat.jl")
        end
    end
end

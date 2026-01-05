"""
    FVMProblem{T,DF,DP,RF,RP,L,R,IC,FT}

Definition of an `FVMProblem`.

# Fields

  - `geometry::FVMGeometry{T}`: The geometry of the problem.
  - `boundary_conditions::BoundaryConditions{L, R}`: The boundary conditions.
  - `diffusion_function::DF`: The diffusion function, of the form `(u, x, t, p) -> Number`, where `p = diffusion_parameters`.
  - `diffusion_parameters::DP`: The parameters for the diffusion function.
  - `reaction_function::RF`: The reaction function, of the form `(u, x, t, p) -> Number`, where `p = reaction_parameters`.
  - `reaction_parameters::RP`: The parameters for the reaction function.
  - `initial_condition::IC`: The initial condition, with `initial_condition[i]` corresponding to the value at `geometry.mesh_points[i]` and `t = initial_time`.
  - `initial_time::FT`: The initial time.
  - `final_time::FT`: The final time.

# Constructors

You can use the default constructor, but we also provide the constructor

    FVMProblem(;
        geometry, 
        boundary_conditions,
        diffusion_function,
        diffusion_parameters = nothing,
        reaction_function = Returns(0.0),
        reaction_parameters = nothing,
        initial_condition,
        initial_time = 0.0,
        final_time)

which provides some default values. Moreover, instead of providing `geometry` and `boundary_conditions`, you can use

    FVMProblem(mesh_points, lhs, rhs; kwargs...)

which will construct `geometry = FVMGeometry(mesh_points)` and `boundary_conditions = BoundaryConditions(lhs, rhs)`.
The `kwargs...` are as above, except without `geometry` and `boundary_conditions` of course.

To solve the `FVMProblem`, just use `solve` as you would in DifferentialEquations.jl. For example,

    sol = solve(prob, Tsit5(), saveat=0.1)

solves the problem with the `Tsit5()` algorithm, saving the solution every `0.1` units of time from `initial_time` up to,
and including, `final_time`.
"""
Base.@kwdef struct FVMProblem{T, DF, DP, RF, RP, L, R, IC, FT}
    geometry::FVMGeometry{T}
    boundary_conditions::BoundaryConditions{L, R}
    diffusion_function::DF
    diffusion_parameters::DP = nothing
    reaction_function::RF = Returns(0.0)
    reaction_parameters::RP = nothing
    initial_condition::IC
    initial_time::FT = 0.0
    final_time::FT
end
function FVMProblem(mesh_points, lhs, rhs; kwargs...)
    return FVMProblem(;
        geometry = FVMGeometry(mesh_points),
        boundary_conditions = BoundaryConditions(lhs, rhs),
        kwargs...
    )
end

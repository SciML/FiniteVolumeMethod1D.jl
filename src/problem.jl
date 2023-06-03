"""
    FVMProblem{T,DF,DP,Dθ,RF,RP,Rθ,IC,FT}

Definition of an `FVMProblem`.

# Fields 
- `geometry::FVMGeometry{T}`: The geometry of the problem.
- `diffusion_function::DF`: The diffusion function.
- `diffusion_parameters::DP`: The parameters for the diffusion function.
- `diffusion_theta::Dθ`: The parameters for the diffusion function that are to be estimated.
- `reaction_function::RF`: The reaction function.
- `reaction_parameters::RP`: The parameters for the reaction function.
- `reaction_theta::Rθ`: The parameters for the reaction function that are to be estimated.
- `initial_condition::IC`: The initial condition.
- `initial_time::FT`: The initial time.
- `final_time::FT`: The final time.

# Constructors

You can use the default constructor, but we also provide the constructor 

    FVMProblem(;
        geometry, 
        diffusion_function,
        diffusion_parameters = nothing,
        diffusion_theta = nothing,
        reaction_function = Returns(0.0),
        reaction_parameters = nothing,
        reaction_theta = nothing,
        initial_condition,
        initial_time = 0.0,
        final_time)

which provides some default values. Moreover, instead of providing `geometry`, you can use 

    FVMProblem(mesh_points; kwargs...)

which will construct `geometry = FVMGeometry(mesh_points)`. The `kwargs...` are as above, except 
without `geometry` of course.

To solve the `FVMProblem`, just use `solve` as you would in DifferentialEquations.jl. For example, 

    sol = solve(prob, Tsit5(), saveat=0.1)
"""
Base.@kwdef struct FVMProblem{T,DF,DP,Dθ,RF,RP,Rθ,IC,FT}
    geometry::FVMGeometry{T}
    diffusion_function::DF
    diffusion_parameters::DP = nothing
    diffusion_theta::Dθ = nothing
    reaction_function::RF = Returns(0.0)
    reaction_parameters::RP = nothing
    reaction_theta::Rθ = nothing
    initial_condition::IC
    initial_time::FT = 0.0
    final_time::FT
end
FVMProblem(mesh_points; kwargs...) = FVMProblem(; geometry=FVMGeometry(mesh_points), kwargs...)
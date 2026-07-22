"""
    FVMProblem(;
        geometry, boundary_conditions, diffusion_function, initial_condition,
        final_time, diffusion_parameters = nothing, reaction_function = Returns(0.0),
        reaction_parameters = nothing, initial_time = 0.0
    )
    FVMProblem(mesh_points, lhs, rhs; kwargs...)

Defines a one-dimensional diffusion-reaction finite-volume problem.

# Arguments

- `geometry`: An [`FVMGeometry`](@ref) describing the spatial mesh.
- `boundary_conditions`: A [`BoundaryConditions`](@ref) instance for the mesh endpoints.
- `diffusion_function`: Function called as `(u, x, t, p)` to evaluate diffusion.
- `initial_condition`: Values at the mesh points at `initial_time`.
- `final_time`: Final integration time.

# Keywords

- `diffusion_parameters = nothing`: Parameters passed to `diffusion_function`.
- `reaction_function = Returns(0.0)`: Function called as `(u, x, t, p)` for the reaction.
- `reaction_parameters = nothing`: Parameters passed to `reaction_function`.
- `initial_time = 0.0`: Initial integration time.

# Fields

- `geometry::FVMGeometry{T}`: Spatial mesh geometry.
- `boundary_conditions::BoundaryConditions{L, R}`: Endpoint boundary conditions.
- `diffusion_function::DF`: Diffusion function called as `(u, x, t, p)`.
- `diffusion_parameters::DP`: Parameters passed to `diffusion_function`.
- `reaction_function::RF`: Reaction function called as `(u, x, t, p)`.
- `reaction_parameters::RP`: Parameters passed to `reaction_function`.
- `initial_condition::IC`: State values at `initial_time`.
- `initial_time::FT`: Initial integration time.
- `final_time::FT`: Final integration time.

# Example

```julia
mesh_points = range(0.0, 1.0; length = 11)
problem = FVMProblem(
    mesh_points,
    Dirichlet(0.0),
    Dirichlet(1.0);
    diffusion_function = (u, x, t, p) -> 1.0,
    initial_condition = collect(mesh_points),
    final_time = 0.1,
)
```
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

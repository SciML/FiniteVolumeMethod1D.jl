using ..FiniteVolumeMethod1D

mesh_points = sort(rand(100))
geo = FVMGeometry(mesh_points)
diffusion_function = (u, p) -> u * p[1] + 2.0
diffusion_parameters = (1.0,)
reaction_function = (u, p) -> u * p[1] + 2.0 + p[2]^2
reaction_parameters = (1.0, 2.0)
initial_condition = rand(100)
final_time = 2.
lhs = Dirichlet((u, p) -> 0.39u + p, 0.5)
rhs = Neumann((u, p) -> 5.8)
boundary_conditions = BoundaryConditions(lhs, rhs)
prob = FVMProblem(;
    geometry=geo,
    boundary_conditions=boundary_conditions,
    diffusion_function=diffusion_function,
    diffusion_parameters=diffusion_parameters,
    reaction_function=reaction_function,
    reaction_parameters=reaction_parameters,
    initial_condition=initial_condition,
    final_time=final_time)
@test prob.geometry == geo
@test prob.boundary_conditions == boundary_conditions
@test prob.boundary_conditions.lhs == lhs
@test prob.boundary_conditions.rhs == rhs
@test prob.diffusion_function == diffusion_function
@test prob.diffusion_parameters == diffusion_parameters
@test prob.reaction_function == reaction_function
@test prob.reaction_parameters == reaction_parameters
@test prob.initial_condition == initial_condition
@test prob.final_time == final_time
@test prob.initial_time == 0.0

prob = FVMProblem(mesh_points, lhs, rhs;
    diffusion_function=diffusion_function,
    diffusion_parameters=diffusion_parameters,
    reaction_function=reaction_function,
    reaction_parameters=reaction_parameters,
    initial_condition=initial_condition,
    final_time=final_time,
    initial_time=0.2)
@test prob.geometry.mesh_points == geo.mesh_points 
@test prob.geometry.spacings == geo.spacings
@test prob.geometry.volumes == geo.volumes
@test prob.boundary_conditions == boundary_conditions
@test prob.boundary_conditions.lhs == lhs
@test prob.boundary_conditions.rhs == rhs
@test prob.diffusion_function == diffusion_function
@test prob.diffusion_parameters == diffusion_parameters
@test prob.reaction_function == reaction_function
@test prob.reaction_parameters == reaction_parameters
@test prob.initial_condition == initial_condition
@test prob.final_time == final_time
@test prob.initial_time == 0.2
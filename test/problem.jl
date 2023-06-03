using ..FiniteVolumeMethod1D

mesh_points = sort(rand(100))
geo = FVMGeometry(mesh_points)
diffusion_function = (u, p) -> u * p[1] + 2.0
diffusion_parameters = (1.0,)
diffusion_theta = [0.1, 0.5]
reaction_function = (u, p) -> u * p[1] + 2.0 + p[2]^2
reaction_parameters = (1.0, 2.0)
initial_condition = rand(100)
final_time = 2.8771
prob = FVMProblem(;
    geometry=geo,
    diffusion_function=diffusion_function,
    diffusion_parameters=diffusion_parameters,
    diffusion_theta=diffusion_theta,
    reaction_function=reaction_function,
    reaction_parameters=reaction_parameters,
    initial_condition=initial_condition,
    final_time=final_time)
@test prob.geometry == geo
@test prob.diffusion_function == diffusion_function
@test prob.diffusion_parameters == diffusion_parameters
@test prob.reaction_function == reaction_function
@test prob.reaction_parameters == reaction_parameters
@test prob.initial_condition == initial_condition
@test prob.final_time == final_time
@test prob.initial_time == 0.0
@test prob.diffusion_theta == diffusion_theta
@test prob.reaction_theta === nothing
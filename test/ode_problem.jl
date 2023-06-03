using ..FiniteVolumeMethod1D
using SciMLBase
using LinearAlgebra
FVM = FiniteVolumeMethod1D

mesh_points = sort(rand(100))
geo = FVMGeometry(mesh_points)
diffusion_function = (u, p) -> u * p[1] + 2.0
diffusion_parameters = (1.0,)
reaction_function = (u, p) -> u * p[1] + 2.0 + p[2]^2
reaction_parameters = (1.0, 2.0)
initial_condition = rand(100)
initial_time = 0.3
final_time = 2.8771
lhs = Dirichlet(0.5)
rhs = Robin((u, p) -> (0.39u + p, 0.3), 0.5)
prob = FVMProblem(;
    geometry=geo,
    boundary_conditions = BoundaryConditions(lhs, rhs),
    diffusion_function=diffusion_function,
    diffusion_parameters=diffusion_parameters,
    reaction_function=reaction_function,
    reaction_parameters=reaction_parameters,
    initial_condition=initial_condition,
    final_time=final_time,
    initial_time=initial_time)
J = FVM.jacobian_sparsity(prob)
@test J == FVM.jacobian_sparsity(mesh_points)
_J = zeros(length(mesh_points), length(mesh_points))
for i in eachindex(mesh_points)
    if i == 1
        idx = [i, i + 1]
    elseif i == lastindex(mesh_points)
        idx = [i - 1, i]
    else
        idx = [i - 1, i, i + 1]
    end
    _J[i, idx] .= 1
end
@test J == _J
@test Tridiagonal(J) == J
ode_prob = ODEProblem(prob)
@test ode_prob.p == prob
@test ode_prob.f.f.f == FVM.pde_odes!
@test ode_prob.u0 == initial_condition
@test ode_prob.tspan == (prob.initial_time, final_time)
@test ode_prob.f.jac_prototype == J
cb = ode_prob.kwargs.data.callback.discrete_callbacks
lhs_cb, rhs_cb = cb 
@test lhs_cb.condition(0.0,0.0,0.0)
@test !rhs_cb.condition(0.0,0.0,0.0)
@test lhs_cb.affect! == FVM.update_lhs!
@test rhs_cb.affect! == FVM.update_rhs!
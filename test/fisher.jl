
using ..FiniteVolumeMethod1D
using MethodOfLines
using OrdinaryDiffEq
using ModelingToolkit
using DomainSets

# Solve: du/dt = d/dx[D(u)du/dx] + R(u),
#   where   D(u) = 1/(10u) + 50/u^2 + 3/u^3
#           R(u) = βKu(1 - u/K)
#           u(0, x) = x < 0 ? 1 : 1/2
#           ∂ₓu(t, -6π) = 0
#           ∂ₓu(t, 6π) = 0,
#           for -6π ≤ x ≤ 6 and 0 ≤ t \leq 5

@parameters t x
@parameters θ₁ θ₂ θ₃
@parameters β K
@variables u(..)
ic_f(x) = x < 0 ? 1.0 : 1 / 2
@register_symbolic ic_f(x)
Dt = Differential(t)
Dx = Differential(x)
Dxx = Differential(x)^2
eqs = [Dt(u(t, x)) ~
    Dx(
        (inv(θ₁ * u(t, x)) + θ₂ * inv(u(t, x)^2) + θ₃ * inv(u(t, x)^3)) *
        Dx(u(t, x))
    ) +
    β * K * u(t, x) * (1 - u(t, x) / K)]
bcs = [Dx(u(t, 2π)) ~ 0.0,
    Dx(u(t, -2π)) ~ 0.0,
    u(0, x) ~ ic_f(x)]
domains = [t ∈ Interval(0.0, 1.0),
    x ∈ Interval(-2π, 2π)]
@named pdesys = PDESystem(
    eqs,
    bcs,
    domains,
    [t, x],
    [u(t, x)],
    [
        θ₁ => 1.0,
        θ₂ => 50.0,
        θ₃ => 3.0,
        β => 1e-3,
        K => 2.0
    ]
)
discretisation = MOLFiniteDifference([x => 0.1], t)
prob = discretize(pdesys, discretisation)
saveat = 0.1
sol = solve(prob, Tsit5(), saveat=saveat)
solt = sol[t]
solx = sol[x]
solu = sol[u(t, x)]

diffusion_function = (u, p) -> inv(p[1] * u) + p[2] * inv(u^2) + p[3] * inv(u^3)
reaction_function = (u, p) -> p[1] * p[2] * u * (1 - u / p[2])
diffusion_parameters = [1.0, 50.0, 3.0]
reaction_parameters = [1e-3, 2.0]
mesh_points = solx
initial_condition = ic_f.(mesh_points)
final_time = 1.0
lhs = Neumann(0.0)
rhs = Neumann(0.0)
fvm_prob = FVMProblem(
    mesh_points,
    lhs,
    rhs;
    diffusion_function,
    reaction_function,
    diffusion_parameters,
    reaction_parameters,
    initial_condition,
    final_time
)
fvm_sol = solve(fvm_prob, TRBDF2(), saveat=solt)
fvm_solu = reduce(hcat, fvm_sol.u)'
@test solu ≈ fvm_solu rtol = 1e-2
using FiniteVolumeMethod1D
using LinearAlgebra
using LinearSolve
using OrdinaryDiffEq
using OrdinaryDiffEqSDIRK

# p. 145 Constanda
# uₜ = uₓₓ + π²exp[-24π²t]sin(5πx), 0 < x < 1
# u(0, t) = 0
# u(1, t) = 0
# u(x, 0) = 3sin(4πx)

function exact_u(x, t)
    return 3exp(-16π^2 * t) * sin(4π * x) + exp(-25π^2 * t) * (exp(π^2 * t) - 1) * sin(5π * x)
end

reaction_function = (u, x, t, p) -> p[1] * exp(-p[2] * t) * sin(p[3] * x)
reaction_parameters = (π^2, 24π^2, 5π)
diffusion_function = (u, x, t, p) -> one(u)
lhs = Dirichlet(0.0)
rhs = Dirichlet(0.0)
mesh_points = LinRange(0, 1, 500)
ic_ff = x -> 3sin(4π * x)
initial_condition = ic_ff.(mesh_points)
final_time = 0.05
prob = FVMProblem(
    mesh_points, lhs, rhs;
    diffusion_function,
    reaction_function,
    reaction_parameters,
    initial_condition,
    final_time
)
sol = solve(
    prob, TRBDF2(linsolve = KLUFactorization()),
    saveat = 0.001, abstol = 1.0e-10, reltol = 1.0e-8
)
exact_sol = [exact_u.(mesh_points, sol.t[i]) for i in eachindex(sol.t)]
@test reduce(hcat, sol.u) ≈ reduce(hcat, exact_sol) rtol = 1.0e-1

using FiniteVolumeMethod1D
using OrdinaryDiffEq
using OrdinaryDiffEqSDIRK
using SpecialFunctions
FVM = FiniteVolumeMethod1D

D₀ = 2.1e-3
u₀ = 8.0e-1
Q = 1.0
m = 2
r₀ = Q * gamma(1 / m + 3 / 2) / (sqrt(π) * u₀ * gamma(1 / m + 1))
t₀ = r₀^2 * m / (2D₀ * (m + 2))
λ = t -> (t / t₀)^(1 / (2 + m))
exact_solution = (
    x, t,
) -> abs(x) > r₀ * λ(t) ? zero(x) : u₀ / λ(t) * (1 - (x / (r₀ * λ(t)))^2)^(1 / m)
L = 15.0
Q = 1.0
N = 2001
mesh_points = LinRange(-L, L, N)
N′ = (N + 1) ÷ 2
base = mesh_points[N′ + 1] - mesh_points[N′ - 1]
initial_condition = zeros(N)
initial_condition[N′] = 2Q / base # Approximating a delta initial condition with a triangle around 0, making its area equal the mass Q by changing the point at 0
initial_time = 0.0
final_time = 5.0
saveat = [0.1, 0.5, 1.0, 2.0, 3.0, 5.0]
diffusion = (u, x, t, p) -> p.p * (u / p.θ[1])^p.θ[2]
diffusion_p = D₀
diffusion_θ = [u₀, m]
diffusion_parameters = (p = diffusion_p, θ = diffusion_θ)
lhs = Neumann(0.0)
rhs = Neumann(0.0)
prob = FVMProblem(
    mesh_points, lhs, rhs;
    initial_time,
    final_time,
    initial_condition,
    diffusion_function = diffusion,
    diffusion_parameters,
    reaction_function = (u, x, t, p) -> zero(u)
)
sol = solve(prob, TRBDF2(); saveat = saveat)
exact = [exact_solution.(mesh_points, sol.t[i]) for i in eachindex(sol.t)]
@test exact ≈ sol.u rtol = 1.0e-1

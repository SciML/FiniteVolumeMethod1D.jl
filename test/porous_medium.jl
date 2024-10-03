using ..FiniteVolumeMethod1D
using OrdinaryDiffEq
using CairoMakie
using ReferenceTests
using SpecialFunctions
FVM = FiniteVolumeMethod1D

D₀ = 2.1e-3
u₀ = 8e-1
Q = 1.0
m = 2
r₀ = Q * gamma(1 / m + 3 / 2) / (sqrt(π) * u₀ * gamma(1 / m + 1))
t₀ = r₀^2 * m / (2D₀ * (m + 2))
λ = t -> (t / t₀)^(1 / (2 + m))
exact_solution = (x, t) -> abs(x) > r₀ * λ(t) ? zero(x) : u₀ / λ(t) * (1 - (x / (r₀ * λ(t)))^2)^(1 / m)
L = 15.0
Q = 1.0
N = 2001
mesh_points = LinRange(-L, L, N)
N′ = (N + 1) ÷ 2
base = mesh_points[N′+1] - mesh_points[N′-1]
initial_condition = zeros(N)
initial_condition[N′] = 2Q / base # Approximating a delta initial condition with a triangle around 0, making its area equal the mass Q by changing the point at 0
initial_time = 0.0
final_time = 5.0
saveat = [0.1, 0.5, 1.0, 2.0, 3.0, 5.0]
diffusion = (u, x, t, p) -> p.p * (u / p.θ[1])^p.θ[2]
diffusion_p = D₀
diffusion_θ = [u₀, m]
diffusion_parameters = (p=diffusion_p, θ=diffusion_θ)
lhs = Neumann(0.0)
rhs = Neumann(0.0)
prob = FVMProblem(mesh_points, lhs, rhs;
    initial_time,
    final_time,
    initial_condition,
    diffusion_function=diffusion,
    diffusion_parameters,
    reaction_function=(u, x, t, p) -> zero(u)
)
sol = solve(prob, TRBDF2(); saveat=saveat)
exact = [exact_solution.(mesh_points, sol.t[i]) for i in eachindex(sol)]
@test exact ≈ sol.u rtol = 1e-1

fig = Figure(size=(2150, 460), fontsize=34)
ax = Axis(fig[1, 1], xlabel=L"x", ylabel=L"u(x)", title=L"(a):$ $ Numerical solution", titlealign=:left, width=600, height=300)
colors = [:red, :black, :blue, :darkgreen, :magenta, :orange]
[lines!(ax, mesh_points, sol.u[i], color=colors[i]) for i in eachindex(sol)]
ylims!(ax, -1e-6, 5)
xlims!(ax, -1, 1)
ax = Axis(fig[1, 2], xlabel=L"x", ylabel=L"u(x)", title=L"(b):$ $ Exact solution", titlealign=:left, width=600, height=300)
[lines!(ax, mesh_points, exact[i], color=colors[i]) for i in eachindex(sol)]
ylims!(ax, -1e-6, 5)
xlims!(ax, -1, 1)
ax = Axis(fig[1, 3], xlabel=L"x", ylabel=L"u(x)", title=L"(c):$ $ Error", titlealign=:left, width=600, height=300)
[lines!(ax, mesh_points, abs.(exact[i] .- sol.u[i]), color=colors[i]) for i in eachindex(sol)]
ylims!(ax, -1e-6, 0.5)
xlims!(ax, -1, 1)
fig_path = normpath(@__DIR__, "..", "test", "figures")
@test_reference joinpath(fig_path, "porous_comparison.png") fig

let t_range = LinRange(0.0, final_time, 250)
    fig = Figure(fontsize=33)
    ax = Axis3(fig[1, 1], xlabel=L"x", ylabel=L"t", zlabel=L"z")
    sol_u = [sol(t) for t in t_range]
    in_int_range = findall(x -> abs(x) < 1, mesh_points) # Makie won't clip the plot to (-1, 1)
    surface!(ax, mesh_points[in_int_range], t_range, reduce(hcat, sol_u)[in_int_range, :], colormap=:viridis)
    xlims!(ax, -1, 1)
    fig_path = normpath(@__DIR__, "..", "docs", "src", "figures")
    @test_reference joinpath(fig_path, "porous_surface.png") fig
end
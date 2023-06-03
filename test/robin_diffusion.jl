using ..FiniteVolumeMethod1D
using LinearAlgebra
using LinearSolve
using OrdinaryDiffEq
using NonlinearSolve
using CairoMakie 
using ReferenceTests

# uₜ = (1/25)uₓₓ
# u(0, t) = 0
# uₓ(3, t) = -u(3, t)/2
# u(x, 0) = 100(1 - x/3)
# u(x, t) ≈ 47.0449exp[-0.0210t]sin(0.7249x) + ⋯
# See http://ramanujan.math.trinity.edu/rdaileda/teach/s12/m3357/lectures/lecture_2_28_short.pdf

# The eigenvalue problem to solve is tan(μₙL) = -μₙ/κ, where L = 3 and κ = 1/2
function compute_μₙ(n)
    interval = ((2n - 1)π / 6+0.001, n * π / 3-0.001)
    f = (μₙ, _) -> tan(3μₙ) + 2μₙ
    prob = IntervalNonlinearProblem(f, interval)
    sol = solve(prob, Ridder(), reltol=1e-9)
    return sol.u
end
function exact_u(x, t, μ)
    u = 0.0
    for (n, μₙ) in enumerate(μ)
        u += 200 * (3μₙ - sin(3μₙ)) * exp(-μₙ^2 * t/25) * sin(μₙ*x) / (3μₙ^2 * (3 + 2cos(3μₙ)^2))
    end
    return u 
end

c = 1 / 5
mesh_points = LinRange(0, 3, 2500)
lhs = Dirichlet(0.0)
rhs = Robin((u, t, p) -> (p[1] * u, p[2]), (1 / 2, 1.0))
diffusion_function = (u, x, t, p) -> p^2
diffusion_parameters = c
ic = x -> 100(1 - x / 3)
initial_condition = ic.(mesh_points)
prob = FVMProblem(mesh_points, lhs, rhs;
    diffusion_function,
    diffusion_parameters,
    final_time=3.0,
    initial_condition
)
sol = solve(prob, TRBDF2(linsolve=KLUFactorization()), saveat = 0.01)

μ = compute_μₙ.(1:100)
exact_sol = [exact_u.(mesh_points, sol.t[i], Ref(μ)) for i in eachindex(sol)]
@test reduce(hcat, sol.u) ≈ reduce(hcat, exact_sol) rtol = 1e-2

let t_range = LinRange(0.0, final_time, 250)
    fig = Figure(fontsize=33)
    ax = Axis3(fig[1, 1], xlabel=L"x", ylabel=L"t", zlabel=L"z", azimuth = 0.8)
    sol_u = [sol(t) for t in t_range]
    surface!(ax, mesh_points, t_range, reduce(hcat, sol_u), colormap=:viridis)
    fig_path = normpath(@__DIR__, "..", "docs", "src", "figures")
    @test_reference joinpath(fig_path, "robin_diffusion_surface.png") fig
end
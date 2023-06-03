using ..FiniteVolumeMethod1D
using OrdinaryDiffEq
using ReferenceTests
using CairoMakie
using LinearSolve

function exact2(x, t)
    u = 1 / 2
    for n in 1:20
        u += ((-1)^n - 1) * 2 * cos(n * π * x) * exp(-n^2 * π^2 * t) / (n^2 * π^2)
    end
    return u
end

mesh_points = LinRange(0, 1, 500)
lhs = Neumann(0.0)
rhs = Neumann(0.0)

diffusion_function = (u, x, t, p) -> one(u)
initial_condition = collect(mesh_points)
final_time = 0.25
prob = FVMProblem(mesh_points, lhs, rhs;
    diffusion_function,
    initial_condition,
    final_time)

sol = solve(prob, TRBDF2(linsolve=KLUFactorization()))
exact_sol = [exact2.(mesh_points, sol.t[i]) for i in eachindex(sol)]
@test reduce(hcat, sol.u) ≈ reduce(hcat, exact_sol) rtol = 1e-2

let t_range = LinRange(0.0, final_time, 250)
    fig = Figure(fontsize=33)
    ax = Axis(fig[1, 1], xlabel=L"x", ylabel=L"t")
    sol_u = [sol(t) for t in t_range]
    contourf!(ax, mesh_points, t_range, reduce(hcat, sol_u), colormap=:viridis)
    fig_path = normpath(@__DIR__, "..", "docs", "src", "figures")
    @test_reference joinpath(fig_path, "heat_contour.png") fig
    fig
end
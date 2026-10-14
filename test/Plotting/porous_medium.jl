using CairoMakie
using ReferenceTests
include(joinpath(@__DIR__, "..", "porous_medium.jl"))

fig = Figure(size = (2150, 460), fontsize = 34)
ax = Axis(
    fig[1, 1], xlabel = L"x", ylabel = L"u(x)", title = L"(a):$ $ Numerical solution",
    titlealign = :left, width = 600, height = 300
)
colors = [:red, :black, :blue, :darkgreen, :magenta, :orange]
[lines!(ax, mesh_points, sol.u[i], color = colors[i]) for i in eachindex(sol.t)]
ylims!(ax, -1.0e-6, 5)
xlims!(ax, -1, 1)
ax = Axis(
    fig[1, 2], xlabel = L"x", ylabel = L"u(x)", title = L"(b):$ $ Exact solution",
    titlealign = :left, width = 600, height = 300
)
[lines!(ax, mesh_points, exact[i], color = colors[i]) for i in eachindex(sol.t)]
ylims!(ax, -1.0e-6, 5)
xlims!(ax, -1, 1)
ax = Axis(
    fig[1, 3], xlabel = L"x", ylabel = L"u(x)", title = L"(c):$ $ Error",
    titlealign = :left, width = 600, height = 300
)
[
    lines!(ax, mesh_points, abs.(exact[i] .- sol.u[i]), color = colors[i])
        for i in eachindex(sol.t)
]
ylims!(ax, -1.0e-6, 0.5)
xlims!(ax, -1, 1)
fig_path = normpath(@__DIR__, "..", "figures")
@test_reference joinpath(fig_path, "porous_comparison.png") fig

let t_range = LinRange(0.0, final_time, 250)
    fig = Figure(fontsize = 33)
    ax = Axis3(fig[1, 1], xlabel = L"x", ylabel = L"t", zlabel = L"z")
    sol_u = [sol(t) for t in t_range]
    in_int_range = findall(x -> abs(x) < 1, mesh_points) # Makie won't clip the plot to (-1, 1)
    surface!(
        ax, mesh_points[in_int_range], t_range,
        reduce(hcat, sol_u)[in_int_range, :], colormap = :viridis
    )
    xlims!(ax, -1, 1)
    fig_path = normpath(@__DIR__, "..", "..", "docs", "src", "figures")
    @test_reference joinpath(fig_path, "porous_surface.png") fig
end

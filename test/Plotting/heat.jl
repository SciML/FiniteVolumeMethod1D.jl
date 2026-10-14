using CairoMakie
using ReferenceTests
include(joinpath(@__DIR__, "..", "heat.jl"))

let t_range = LinRange(0.0, final_time, 250)
    fig = Figure(fontsize = 33)
    ax = Axis(fig[1, 1], xlabel = L"x", ylabel = L"t")
    sol_u = [sol(t) for t in t_range]
    contourf!(ax, mesh_points, t_range, reduce(hcat, sol_u), colormap = :viridis)
    fig_path = normpath(@__DIR__, "..", "..", "docs", "src", "figures")
    @test_reference joinpath(fig_path, "heat_contour.png") fig
    fig
end

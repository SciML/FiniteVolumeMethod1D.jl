using CairoMakie
using ReferenceTests
include(joinpath(@__DIR__, "..", "robin_diffusion.jl"))

let t_range = LinRange(0.0, final_time, 250)
    fig = Figure(fontsize = 33)
    ax = Axis3(fig[1, 1], xlabel = L"x", ylabel = L"t", zlabel = L"z", azimuth = 0.8)
    sol_u = [sol(t) for t in t_range]
    surface!(ax, mesh_points, t_range, reduce(hcat, sol_u), colormap = :viridis)
    fig_path = normpath(@__DIR__, "..", "..", "docs", "src", "figures")
    @test_reference joinpath(fig_path, "robin_diffusion_surface.png") fig
end

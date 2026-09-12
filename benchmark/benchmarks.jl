using FiniteVolumeMethod1D, BenchmarkTools
using OrdinaryDiffEqSDIRK, LinearSolve, StableRNGs
using LinearAlgebra

const SUITE = BenchmarkGroup()
const rng = StableRNG(123)

# =============================================================================
# Geometry and boundary conditions
# =============================================================================

mesh_points = LinRange(0, 1, 500)
ic_ff = x -> 3sin(4π * x)
initial_condition = ic_ff.(mesh_points)
final_time = 0.05

diffusion_function = (u, x, t, p) -> one(u)
reaction_function = (u, x, t, p) -> p[1] * exp(-p[2] * t) * sin(p[3] * x)
reaction_parameters = (π^2, 24π^2, 5π)

lhs = Dirichlet(0.0)
rhs = Dirichlet(0.0)

SUITE["construct"] = BenchmarkGroup()

SUITE["construct"]["geometry"] = @benchmarkable FVMGeometry($mesh_points)
SUITE["construct"]["dirichlet"] = @benchmarkable Dirichlet(0.0)
SUITE["construct"]["neumann"] = @benchmarkable Neumann(0.5)
SUITE["construct"]["bcs"] = @benchmarkable BoundaryConditions($lhs, $rhs)
SUITE["construct"]["problem"] = @benchmarkable FVMProblem(
    $mesh_points, $lhs, $rhs; diffusion_function = $diffusion_function,
    reaction_function = $reaction_function,
    reaction_parameters = $reaction_parameters,
    initial_condition = $initial_condition, final_time = $final_time
)

# =============================================================================
# PDE solve — diffusion-reaction with Dirichlet boundaries
# =============================================================================

prob = FVMProblem(
    mesh_points, lhs, rhs; diffusion_function = diffusion_function,
    reaction_function = reaction_function,
    reaction_parameters = reaction_parameters,
    initial_condition = initial_condition, final_time = final_time
)

SUITE["solve"] = BenchmarkGroup()

SUITE["solve"]["trbdf2_klu"] = @benchmarkable solve(
    $prob, TRBDF2(; linsolve = KLUFactorization())
)

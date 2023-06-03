using ..FiniteVolumeMethod1D
using LinearAlgebra
using LinearSolve
using OrdinaryDiffEq

# p. 145 Constanda 
# uₜ = uₓₓ + π²exp[-24π²t]sin(5πx), 0 < x < 1 
# u(0, t) = 0
# u(1, t) = 0
# u(x, 0) = 3sin(4πx)

function exact_u(x, t)
    return 3exp(-16π^2*t)*sin(4π*x) + exp(-25π^2*t)*(exp(π^2*t)-1)*sin(5π*x)
end

reaction_function = (u, p) -> p[1] * exp(-p[2])
diffusion_function = (u, p) -> one(u)
lhs = Diffusion(0.0)
rhs = Diffusion(0.0)
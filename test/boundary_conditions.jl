using ..FiniteVolumeMethod1D
FVM = FiniteVolumeMethod1D

f = (u, t, p) -> 1.01u^p[1] + p[2] + t
p = (1.0, 2.0)
dirc = Dirichlet(f, p)
@test FVM.is_dirichlet(dirc)
@test !FVM.is_neumann(dirc)
@test dirc(0.881, 0.5) == f(0.881, 0.5, p)
@inferred dirc(0.881, 0.5)
@test Dirichlet(f) == Dirichlet(f, nothing)
dirc = Dirichlet(0.5)
@test dirc(0.5, 0.2) == 0.5
@inferred dirc(0.5, 0.2)
dirc = Dirichlet(0)
@test dirc(0, 0.2) == 0 && dirc(0, 0.2) isa Int
@test dirc(0.0, 0.2) == 0.0 && dirc(0.0, 0.2) isa Float64
@inferred dirc(0, 0.2)
@inferred dirc(0.0, 0.2)

neum = Neumann(f, p)
@test !FVM.is_dirichlet(neum)
@test FVM.is_neumann(neum)
@test neum(0.5, 0.25) == f(0.5, 0.25, p)
@test neum(0.881, 0.192) == f(0.881, 0.192, p)
@inferred neum(0.881, 0.192)
@test Neumann(f) == Neumann(f, nothing)
neum = Neumann(0.5)
@test neum(0.5, 0.2) == 0.5
@inferred neum(0.5, 0.2)
neum = Neumann(0)
@test neum(0, 0.2) == 0 && neum(0, 0.2) isa Int
@test neum(0.0, 0.2) == 0.0 && neum(0.0, 0.2) isa Float64
@inferred neum(0, 0.2)

bc = BoundaryConditions(dirc, neum)
@test bc.lhs == dirc
@test bc.rhs == neum

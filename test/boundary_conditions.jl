using ..FiniteVolumeMethod1D
FVM = FiniteVolumeMethod1D

f = (u, p) -> 1.01u^p[1] + p[2]
p = (1.0, 2.0)
dirc = Dirichlet(f, p)
@test FVM.is_dirichlet(dirc)
@test !FVM.is_neumann(dirc)
@test !FVM.is_robin(dirc)
@test_throws ArgumentError FVM.get_ab(dirc, 0.5)
@test dirc(0.881) == f(0.881, p)
@inferred dirc(0.881) 
@test Dirichlet(f) == Dirichlet(f, nothing)
dirc = Dirichlet(0.5)
@test dirc(0.5) == 0.5
@inferred dirc(0.5)
dirc = Dirichlet(0)
@test dirc(0) == 0 && dirc(0) isa Int
@test dirc(0.0) == 0.0 && dirc(0.0) isa Float64
@inferred dirc(0)
@inferred dirc(0.0)


neum = Neumann(f, p)
@test !FVM.is_dirichlet(neum)
@test FVM.is_neumann(neum)
@test !FVM.is_robin(neum)
@test FVM.get_ab(neum, 0.5) == (-f(0.5, p), one(0.5))
@test neum(0.881) == f(0.881, p)
@inferred neum(0.881)
@test Neumann(f) == Neumann(f, nothing)
neum = Neumann(0.5)
@test neum(0.5) == 0.5
@inferred neum(0.5)
neum = Neumann(0)
@test neum(0) == 0 && neum(0) isa Int
@test neum(0.0) == 0.0 && neum(0.0) isa Float64
@inferred neum(0)

f = (u, p) -> (u, p[1] + p[2])
robin = Robin(f, p)
@test !FVM.is_dirichlet(robin)
@test !FVM.is_neumann(robin)
@test FVM.is_robin(robin)
@test FVM.get_ab(robin, 0.5) == (0.5, 3.0)
@test robin(0.881) == (0.881, 3.0)
@inferred robin(0.881)
@test Robin(f) == Robin(f, nothing)
robin = Robin((0.5, 3.0))
@test robin(0.5) == (0.5, 3.0)
@inferred robin(0.5)
robin = Robin(0, 3)
@test robin(0) == (0, 3) && robin(0) isa Tuple{Int, Int}
@test robin(0.3) == (0.0, 3.0) && robin(0.3) isa Tuple{Float64, Float64}
@inferred robin(0)
@inferred robin(0.0)

bc = BoundaryConditions(dirc, neum)
@test bc.lhs == dirc 
@test bc.rhs == neum
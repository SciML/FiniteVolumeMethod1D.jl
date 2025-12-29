using ExplicitImports
using FiniteVolumeMethod1D
using Test

@test check_no_implicit_imports(FiniteVolumeMethod1D) === nothing
@test check_no_stale_explicit_imports(FiniteVolumeMethod1D) === nothing

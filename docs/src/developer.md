```@meta
CurrentModule = FiniteVolumeMethod1D
```

# Developer API

This page documents the extension points used by the finite-volume
boundary-condition kernel. They are versioned developer interfaces, not the
recommended user-facing API. Application code should use [`Dirichlet`](@ref),
[`Neumann`](@ref), and [`BoundaryConditions`](@ref).

```@docs
FiniteVolumeMethod1D.AbstractBoundaryCondition
FiniteVolumeMethod1D.is_dirichlet
FiniteVolumeMethod1D.is_neumann
```

```@meta
CurrentModule = FiniteVolumeMethod1D
```

# API

FiniteVolumeMethod1D defines the mesh, endpoint boundary conditions, and problem
description.
The package also reexports `solve` from CommonSolve so an `FVMProblem` can be solved with a
CommonSolve-compatible solver package.

```@docs
FVMGeometry
BoundaryConditions
Dirichlet
Neumann
FVMProblem
solve
```

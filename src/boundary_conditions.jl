abstract type AbstractBoundaryCondition{F, P} end
(bc::AbstractBoundaryCondition{F, P})(u, t) where {F, P} = bc.f(u, t, bc.p)

"""
    Dirichlet(f, p = nothing)
    Dirichlet(; f, p = nothing)
    Dirichlet(value::Number)

A Dirichlet boundary condition for an [`FVMProblem`](@ref).

`f` must accept `(u, t, p)` and return the prescribed value at a boundary. Passing a
number constructs a constant boundary condition. `p` stores optional parameters passed to
`f`.

# Fields

- `f`: Function called as `f(u, t, p)`.
- `p`: Parameters passed to `f`.

# Example

```julia
left_boundary = Dirichlet(0.0)
right_boundary = Dirichlet((u, t, p) -> p * sin(t), 1.0)
```
"""
Base.@kwdef struct Dirichlet{F, P} <: AbstractBoundaryCondition{F, P}
    f::F
    p::P = nothing
    Dirichlet(f::F, p::P = nothing) where {F, P} = new{F, P}(f, p)
end
Dirichlet(f::Function) = Dirichlet(f, nothing)
Dirichlet(v::Number) =
let v = v
    Dirichlet((u, t, p) -> oftype(u, v))
end

"""
    Neumann(f, p = nothing)
    Neumann(; f, p = nothing)
    Neumann(value::Number)

A Neumann boundary condition for an [`FVMProblem`](@ref).

`f` must accept `(u, t, p)` and return the boundary derivative. Passing a number constructs
a constant derivative condition. `p` stores optional parameters passed to `f`.

# Fields

- `f`: Function called as `f(u, t, p)`.
- `p`: Parameters passed to `f`.

# Example

```julia
left_boundary = Neumann(0.0)
right_boundary = Neumann((u, t, p) -> p * u, -0.5)
```
"""
Base.@kwdef struct Neumann{F, P} <: AbstractBoundaryCondition{F, P}
    f::F
    p::P = nothing
    Neumann(f::F, p::P = nothing) where {F, P} = new{F, P}(f, p)
end
Neumann(v::Number) =
let v = v
    Neumann((u, t, p) -> oftype(u, v))
end

is_dirichlet(::AbstractBoundaryCondition) = false
is_dirichlet(::Dirichlet) = true
is_neumann(::AbstractBoundaryCondition) = false
is_neumann(::Neumann) = true

"""
    BoundaryConditions(lhs, rhs)
    BoundaryConditions(; lhs, rhs)

Stores the left and right boundary conditions of an [`FVMProblem`](@ref).

# Fields

- `lhs::L`: Boundary condition at the first mesh point.
- `rhs::R`: Boundary condition at the last mesh point.

# Example

```julia
boundary_conditions = BoundaryConditions(Dirichlet(0.0), Neumann(0.0))
```

See also [`Dirichlet`](@ref) and [`Neumann`](@ref) for the types of
boundary conditions you can construct.
"""
Base.@kwdef struct BoundaryConditions{L, R}
    lhs::L
    rhs::R
end

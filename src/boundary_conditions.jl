"""
    AbstractBoundaryCondition{F, P}

Developer interface for endpoint boundary conditions used by [`FVMProblem`](@ref).
The package supplies the generic two-argument call `bc(u, t)`, which evaluates the
stored three-argument function `bc.f(u, t, bc.p)`.

# Extension Interface

A custom boundary condition should subtype `AbstractBoundaryCondition{F, P}`, store
fields named `f` and `p`, and make `f(u, t, p)` return the boundary value or flux.
Define `is_dirichlet(::YourBoundaryCondition)` to return `true` for a prescribed-value
condition. Leave it at the default `false` for a flux condition, or define
`is_neumann(::YourBoundaryCondition)` when the classification should be explicit.

These hooks are developer API. End users should normally use [`Dirichlet`](@ref),
[`Neumann`](@ref), and [`BoundaryConditions`](@ref).
"""
abstract type AbstractBoundaryCondition{F, P} end
(bc::AbstractBoundaryCondition{F, P})(u, t) where {F, P} = bc.f(u, t, bc.p)

"""
    is_dirichlet(bc::AbstractBoundaryCondition)

Return whether `bc` imposes a prescribed value at the boundary.

Developer subtypes should specialize this predicate to return `true` for their
Dirichlet-like conditions. The default is `false`.
"""
is_dirichlet(::AbstractBoundaryCondition) = false

"""
    is_neumann(bc::AbstractBoundaryCondition)

Return whether `bc` imposes a prescribed boundary flux.

The default is `false`; developer subtypes may specialize this predicate when an
explicit classification is useful. The finite-volume kernel treats a boundary as
flux-based whenever `is_dirichlet(bc)` is `false`.
"""
is_neumann(::AbstractBoundaryCondition) = false

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

# Returns

- `Dirichlet`: A callable boundary condition whose value is imposed at the endpoint.

# Examples

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

# Returns

- `Neumann`: A callable boundary condition whose flux is used at the endpoint.

# Examples

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

is_dirichlet(::Dirichlet) = true
is_neumann(::Neumann) = true

"""
    BoundaryConditions(lhs, rhs)
    BoundaryConditions(; lhs, rhs)

Stores the left and right boundary conditions of an [`FVMProblem`](@ref).

# Fields

- `lhs::L`: Boundary condition at the first mesh point.
- `rhs::R`: Boundary condition at the last mesh point.

# Returns

- `BoundaryConditions`: A pair of endpoint boundary conditions.

# Examples

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

abstract type AbstractBoundaryCondition{F,P} end
(bc::AbstractBoundaryCondition{F,P})(u, t) where {F,P} = bc.f(u, t, bc.p)

@doc raw"""
    Dirichlet{F,P} <: AbstractBoundaryCondition{F,P}

A Dirichlet boundary condition with fields `f` and `p` (default `p = nothing`),
with `f` being a function of the form `f(u, p)` and `p` being
the parameters for `f`. 

A Dirichlet boundary condition takes the form

```math
u(a, t) ↤ f(u(a, t), t, p),
```

where `a` is one of the endpoints. 

# Constructors 

    Dirichlet(f::Function, p = nothing) -> Dirichlet(f, p)
    Dirichlet(; f, p = nothing)         -> Dirichlet(f, p)
    Dirichlet(v::Number)                -> Dirichlet((u, t, p) -> oftype(u, v), nothing)
"""
Base.@kwdef struct Dirichlet{F,P} <: AbstractBoundaryCondition{F,P}
    f::F
    p::P = nothing
    Dirichlet(f::F, p::P=nothing) where {F,P} = new{F,P}(f, p)
end
Dirichlet(f::Function) = Dirichlet(f, nothing)
Dirichlet(v::Number) =
    let v = v
        Dirichlet((u, t, p) -> oftype(u, v))
    end

@doc raw"""
    Neumann{F,P} <: AbstractBoundaryCondition{F,P}

A Neumann boundary condition with fields `f` and `p` (default `p = nothing`),
with `f` being a function of the form `f(u, t, p)` and `p` being
the parameters for `f`.

A Neumann boundary condition takes the form

```math
\dfrac{\partial u}{\partial x}(a, t) = f(u(a, t), t, p),
```

where `a` is one of the endpoints. 

# Constructors 

    Neumann(f::Function, p = nothing) -> Neumann(f, p)
    Neumann(; f, p = nothing)         -> Neumann(f, p)
    Neumann(v::Number)                -> Neumann((u, t, p) -> oftype(u, v), nothing)
"""
Base.@kwdef struct Neumann{F,P} <: AbstractBoundaryCondition{F,P}
    f::F
    p::P = nothing
    Neumann(f::F, p::P=nothing) where {F,P} = new{F,P}(f, p)
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
    BoundaryConditions{L, R}

The boundary conditions for the FVMProblem. 
    
# Fields 
- `lhs::L`: The left-hand side boundary condition.
- `rhs::R`: The right-hand side boundary condition.

See also [`Dirichlet`](@ref) and [`Neumann`](@ref) for the types of 
boundary conditions you can construct.
"""
Base.@kwdef struct BoundaryConditions{L,R}
    lhs::L
    rhs::R
end
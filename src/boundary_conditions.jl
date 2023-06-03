abstract type AbstractBoundaryCondition{F,P} end
(bc::AbstractBoundaryCondition{F,P})(u) where {F,P} = bc.f(u, bc.p)

raw"""
    Dirichlet{F,P} <: AbstractBoundaryCondition{F,P}

A Dirichlet boundary condition with fields `f` and `p` (default `p = nothing`),
with `f` being a function of the form `f(u, p)` and `p` being
the parameters for `f`. 

A Dirichlet boundary condition takes the form

```math
u(a, t) = f(u(a, t), a, t, p),
```

where `a` is one of the endpoints. Note that the `u(a, t)` on the left is the new value, and `u(a, t)` on 
the right is the current value, i.e. this is not an implicit boundary condition.

# Constructors 

    Dirichlet(f::Function, p = nothing) -> Dirichlet(f, p)
    Dirichlet(; f, p = nothing)         -> Dirichlet(f, p)
    Dirichlet(v::Number)                -> Dirichlet((u, p) -> oftype(u, v), nothing)
"""
Base.@kwdef struct Dirichlet{F,P} <: AbstractBoundaryCondition{F,P}
    f::F
    p::P = nothing
end
Dirichlet(f::Function) = Dirichlet(f, nothing)
Dirichlet(v::Number) =
    let v = v
        Dirichlet((u, p) -> oftype(u, v))
    end

raw"""
    Neumann{F,P} <: AbstractBoundaryCondition{F,P}

A Neumann boundary condition with fields `f` and `p` (default `p = nothing`),
with `f` being a function of the form `f(u, p)` and `p` being
the parameters for `f`.

A Neumann boundary condition takes the form

```math
\dfrac{\partial u}{\partial x}(a, t) = f(u(a, t), p),
```

where `a` is one of the endpoints. 

# Constructors 

    Neumann(f::Function, p = nothing) -> Neumann(f, p)
    Neumann(; f, p = nothing)         -> Neumann(f, p)
    Neumann(v::Number)                -> Neumann((u, p) -> oftype(u, v), nothing)
"""
Base.@kwdef struct Neumann{F,P} <: AbstractBoundaryCondition{F,P}
    f::F
    p::P = nothing
end
Neumann(f::Function) = Neumann(f, nothing)
Neumann(v::Number) =
    let v = v
        Neumann((u, p) -> oftype(u, v))
    end

raw"""
    Robin{F,P} <: AbstractBoundaryCondition{F,P}

A Robin boundary condition with fields `f` and `p` (default `p = nothing`),
with `f` being a function of the form `f(u, p)` and `p` being
the parameters for `f`. The function `f` should return a `Tuple` of the form 
`(a₀, b₀)`, where `a₀` and `b₀` are the coefficients in the Robin boundary condition.

A Robin boundary condition takes the form

```math
a₀(u(a, t), p) + b₀\dfrac{\partial u}{\partial x}(a, t) = 0,
```

where `a` is one of the endpoints. 

# Constructors 

    Robin(f::Function, p = nothing)             -> Robin(f, p)
    Robin(; f, p = nothing)                     -> Robin(f, p)
    Robin(a::Number, b::Number)                 -> Robin((u, p) -> (oftype(u, a), oftype(u, b)), nothing)
"""
Base.@kwdef struct Robin{F,P} <: AbstractBoundaryCondition{F,P}
    f::F
    p::P = nothing
end
Robin(f::Function) = Robin(f, nothing)
Robin(a::Number, b::Number) =
    let a = a, b = b
        Robin((u, p) -> (oftype(u, a), oftype(u, b)))
    end

is_dirichlet(::AbstractBoundaryCondition) = false
is_dirichlet(::Dirichlet) = true
is_neumann(::AbstractBoundaryCondition) = false
is_neumann(::Neumann) = true
is_robin(::AbstractBoundaryCondition) = false
is_robin(::Robin) = true

get_ab(bc::Dirichlet, u) = throw(ArgumentError("Cannot get a or b for Dirichlet boundary condition."))
get_ab(bc::Neumann, u) = (-bc(u), one(u))
get_ab(bc::Robin, u) =
    let val = bc(u)
        return (val[1], val[2])
    end

"""
    BoundaryConditions{L<:AbstractBoundaryCondition,R<:AbstractBoundaryCondition}

The boundary conditions for the FVMProblem. 
    
# Fields 
- `lhs::L`: The left-hand side boundary condition.
- `rhs::R`: The right-hand side boundary condition.

See also [`Dirichlet`](@ref), [`Neumann`](@ref), and [`Robin`](@ref) for the types of 
boundary conditions you can construct.
"""
Base.@kwdef struct BoundaryConditions{L<:AbstractBoundaryCondition,R<:AbstractBoundaryCondition}
    lhs::L
    rhs::R
end
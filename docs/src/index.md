```@meta
CurrentModule = FiniteVolumeMethod1D
```

# FiniteVolumeMethod1D

Documentation for [FiniteVolumeMethod1D](https://github.com/SciML/FiniteVolumeMethod1D.jl).

This is a package for solving equations of the form

```math
\dfrac{\partial u(x, t)}{\partial t} = \dfrac{\partial}{\partial x}\left(D\left(u, x, t\right)\dfrac{\partial u(x, t)}{\partial x}\right) + R(u, x, t),
```

using the finite volume method over intervals $a \leq x \leq b$ and $t_0 \leq t \leq t_1$, with support for the following types of boundary conditions (shown at $x = a$, but you can mix boundary condition types, e.g. Neumann at $x=a$ and Robin at $x=b$):

```math
\begin{align*}
\begin{array}{rrcl}
\text{Neumann}: & \dfrac{\partial u(a, t)}{\partial x} & = & a_0\left(u(a, t), t\right), \\[9pt]
\text{Dirichlet}: & u(a, t) & = & a_0\left(u(a, t), t\right),
\end{array}
\end{align*}
```

where the Dirichlet condition has $u(a, t)$ mapping from $a_0(u(a, t), t)$ (i.e., it is not an implicit equation for $u(a, t)$).

More information is given in the sidebar, and the docstrings are below.

If you want a more complete two-dimensional version, please see my other package [FiniteVolumeMethod.jl](https://github.com/SciML/FiniteVolumeMethod.jl).

```@index
```

```@autodocs
Modules = [FiniteVolumeMethod1D]
```

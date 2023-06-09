# FiniteVolumeMethod1D

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://DanielVandH.github.io/FiniteVolumeMethod1D.jl/dev/)
[![Build Status](https://github.com/DanielVandH/FiniteVolumeMethod1D.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/DanielVandH/FiniteVolumeMethod1D.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/DanielVandH/FiniteVolumeMethod1D.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/DanielVandH/FiniteVolumeMethod1D.jl)

This is a package for solving equations of the form

$$
\frac{\partial u}{\partial t} = \frac{\partial}{\partial x}\left(D(u, x, t)\frac{\partial u}{\partial x}\right) + R(u, x, t)
$$

using the finite volume method over intervals $a \leq x \leq b$ and $t_0 \leq t \leq t_1$, with support for the following types of boundary conditions (shown at $x = a$, but you can mix boundary condition types, e.g. Neumann at $x=a$ and Dirichlet at $x=b$):

- `Neumann`: $\dfrac{\partial u(a, t)}{\partial x} = a_0\left(u(a, t), t\right)$.
- `Dirichlet`: $u(a, t) = a_0\left(u(a, t), t\right)$ (this is not an implicit equation for $u(a, t)$, rather $u(a, t)$ is mapped from $a_0\left(u(a, t), a, t\right)$.

For examples on how to use it, please see the docs. If you want a more complete two-dimensional version, please see my other package [FiniteVolumeMethod.jl](https://github.com/DanielVandH/FiniteVolumeMethod.jl).

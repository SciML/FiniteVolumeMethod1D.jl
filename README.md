# FiniteVolumeMethod1D

> Due to unforeseen circumstances, I will no longer be able to actively maintain this package. If you are willing and able to help maintain it, please contact me at ![email](https://github.com/DanielVandH/FiniteVolumeMethod.jl/assets/95613936/0cceed96-a640-4356-a78a-c5cb95944b6c). I will still be around to view Issues and review / accept (hopefully small) PRs in the interim.

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://DanielVandH.github.io/FiniteVolumeMethod1D.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://DanielVandH.github.io/FiniteVolumeMethod1D.jl/dev/)
[![Build Status](https://github.com/DanielVandH/FiniteVolumeMethod1D.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/DanielVandH/FiniteVolumeMethod1D.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/DanielVandH/FiniteVolumeMethod1D.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/DanielVandH/FiniteVolumeMethod1D.jl)
[![DOI](https://zenodo.org/badge/648601161.svg)](https://zenodo.org/badge/latestdoi/648601161)

This is a package for solving equations of the form

$$
\frac{\partial u}{\partial t} = \frac{\partial}{\partial x}\left(D(u, x, t)\frac{\partial u}{\partial x}\right) + R(u, x, t)
$$

using the finite volume method over intervals $a \leq x \leq b$ and $t_0 \leq t \leq t_1$, with support for the following types of boundary conditions (shown at $x = a$, but you can mix boundary condition types, e.g. Neumann at $x=a$ and Dirichlet at $x=b$):

- `Neumann`: $\dfrac{\partial u(a, t)}{\partial x} = a_0\left(u(a, t), t\right)$.
- `Dirichlet`: $u(a, t) = a_0\left(u(a, t), t\right)$ (this is not an implicit equation for $u(a, t)$, rather $u(a, t)$ is mapped from $a_0\left(u(a, t), a, t\right)$.

For examples on how to use it, please see the docs. If you want a more complete two-dimensional version, please see my other package [FiniteVolumeMethod.jl](https://github.com/DanielVandH/FiniteVolumeMethod.jl).

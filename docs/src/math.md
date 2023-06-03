```@meta
CurrentModule = FiniteVolumeMethod1D
```

# Mathematical Details 

In this section, we provide some of the mathematical details for discretising the PDEs we consider. Recall that the problems we consider are

```math
\begin{align*}
\dfrac{\partial u}{\partial t} &= \dfrac{\partial}{\partial x}\left(D(u)\dfrac{\partial u}{\partial x}\right) + R(u), & a \leq x \leq b,\, t_0 < t \leq t_1, \\
a_0\left(u(a, t)\right) + b_0\left(u(a, t)\right)\dfrac{\partial u(a, t)}{\partial x} &= 0, & t_0 < t \leq t_1, \\
a_1\left(u(b, t)\right) + b_1\left(u(b, t)\right)\dfrac{\partial u(b, t)}{\partial x} & =0, & t_0 < t \leq t_1, \\
u(x, 0) &= f(x), & a \leq x \leq b.
\end{align*}
```

This is for the Robin boundary condition form at both ends. We assume that $b_0, b_1 \neq 0$. (We also support functions with arguments $x$ and $t$, e.g. $D(u, x, t)$, but for simplicity we omit the $x$ and $t$ arguments.)

## Interior Discretisation 

Let us start by focusing on the discretisation of the PDE itself. We start by defining some grid $x_1, \ldots, x_n$ for the mesh points, assuming $a = x_1 < x_2 < \cdots < x_n = b$. The control volumes are defined by intervals $[w_i, e_i]$, where

```math
\begin{align*}
w_i &= \begin{cases} x_1 & i=1, \\ \dfrac{1}{2}\left(x_{i-1} + x_i\right) & i=2,\ldots,n, \end{cases} \\
e_i &= \begin{cases} \frac12\left(x_i + x_{i+1}\right) & i=1,\ldots,n-1, \\ x_n & i=n. \end{cases} 
\end{align*}
```

The volumes of these control volumes are defined by $V_i = e_i - w_i$, $i=1,\ldots,n$. Now, integrate the PDE over a control volume $[w_i, e_i]$:

```math
\begin{align*}
\int_{w_i}^{e_i} \dfrac{\partial u}{\partial t}\,\mathrm{d}x &= \int_{w_i}^{e_i} \dfrac{\partial}{\partial x}\left(D(u)\dfrac{\partial u}{\partial x}\right)\,\mathrm{d}x + \int_{w_i}^{e_i} R(u)\,\mathrm{d}x \\
\dfrac{\mathrm d}{\mathrm dt}\int_{w_i}^{e_i} u\,\mathrm{d}x &= D\left(u(e_i, t)\right)\dfrac{\partial u(e_i, t)}{\partial x} - D\left(u(w_i, t)\right)\dfrac{\partial u(w_i, t)}{\partial x} + \int_{w_i}^{e_i} R(u)\,\mathrm{d}x \\
\dfrac{\mathrm d\bar u_i}{\mathrm dt} &= \frac{1}{V_i}\left[D\left(u(e_i, t)\right)\dfrac{\partial u(e_i, t)}{\partial x} - D\left(u(w_i, t)\right)\dfrac{\partial u(w_i, t)}{\partial x}\right] + \bar R_i,
\end{align*}
```

where $\bar u_i = (1/V_i)\int_{w_i}^{e_i} u\,\mathrm{d}x$ and $\bar R_i = (1/V_i)\int_{w_i}^{e_i} R(u)\,\mathrm{d}x$. Letting $u_i = u(x_i, t)$ and $R_i = R(u_i)$, we make the following approximations:

```math
\begin{align*}
\begin{array}{rcll}
\bar u_i &\approx& u_i, & i=1,\ldots,n, \\[6pt]
\bar R_i &\approx& R_i, & i=1,\ldots,n, \\[6pt]
D\left(u(e_i, t)\right) &\approx& \dfrac12\left(D_i + D_{i+1}\right)\quad & i=1,\ldots,n-1, \\[6pt]
D\left(u(w_i, t)\right) &\approx& \dfrac12\left(D_{i-1} + D_i\right)\quad & i=2,\ldots,n, \\[6pt]
\dfrac{\partial u(e_i, t)}{\partial x} &\approx& \dfrac{u_{i+1} - u_i}{h_i} & i=1,\ldots,n-1, \\[6pt]
\dfrac{\partial u(w_i, t)}{\partial x} &\approx& \dfrac{u_i - u_{i-1}}{h_{i-1}} & i=2,\ldots,n, 
\end{array}
\end{align*}
```

where $h_i = x_{i+1} - x_i$, $i=1,\ldots,n-1$. With these approximations, we find:

```math 
\begin{align*}
\frac{\mathrm du_i}{\mathrm dt} &= \frac{1}{V_i}\left[\left(\dfrac{D_i+D_{i+1}}{2}\right)\left(\dfrac{u_{i+1} - u_i}{h_i}\right) - \left(\dfrac{D_{i-1} + D_i}{2}\right)\left(\dfrac{u_i - u_{i-1}}{h_{i-1}}\right)\right] + R_i,
\end{align*}
```

for $i=2,\ldots,n-1$.

## Boundary Discretisation

We still need to handle the equations at $i=1$ and $i=n$. If we rearrange the boundary condition at $x = a$, we obtain

```math
\dfrac{\partial u(a, t)}{\partial x} = -\frac{a_0\left(u(a, t)\right)}{b_0\left(u(a, t)\right)}.
```

Thus,

```math
\dfrac{\mathrm du_1}{\mathrm dt} = \frac{1}{V_1}\left[\left(\dfrac{D_1 + D_2}{2}\right)\left(\dfrac{u_2 - u_1}{h_1}\right) + D(u_1)\frac{a_0(u_1)}{b_0(u_1)}\right] + R_1.
```

Similarly, the boundary condition at $x = b$ gives 

```math 
\dfrac{\partial u(b, t)}{\partial x} = -\dfrac{a_1\left(u(b, t)\right)}{b_1\left(u(b, t)\right)},
```

so

```math
\dfrac{\mathrm du_n}{\mathrm dt} = -\frac{1}{V_n}\left[D(u_n)\frac{a_1(u_n)}{b_1(u_n)} + \left(\dfrac{D_{n-1} + D_n}{2}\right)\left(\dfrac{u_n - u_{n-1}}{h_{n-1}}\right)\right] + R_n.
```

## The Complete Discretisation

Putting all the results together, the complete system of ODEs is

```math
\begin{align*}
\frac{\mathrm du_i}{\mathrm dt} &= \frac{1}{V_i}\left[\left(\dfrac{D_i+D_{i+1}}{2}\right)\left(\dfrac{u_{i+1} - u_i}{h_i}\right) - \left(\dfrac{D_{i-1} + D_i}{2}\right)\left(\dfrac{u_i - u_{i-1}}{h_{i-1}}\right)\right] + R_i,~ i=2,\ldots,n-1, \\[8pt]
\dfrac{\mathrm du_1}{\mathrm dt} &= \frac{1}{V_1}\left[\left(\dfrac{D_1 + D_2}{2}\right)\left(\dfrac{u_2 - u_1}{h_1}\right) + D(u_1)\frac{a_0(u_1)}{b_0(u_1)}\right] + R_1,\\[8pt]
\dfrac{\mathrm du_n}{\mathrm dt}&= \frac{1}{V_n}\left[D(u_n)\frac{a_1(u_n)}{b_1(u_n)} - \left(\dfrac{D_{n-1} + D_n}{2}\right)\left(\dfrac{u_n - u_{n-1}}{h_{n-1}}\right)\right] + R_n.
\end{align*}
```

This system can then be easily solved using methods from DifferentialEquation.jl, treating the system in the form $\boldsymbol u(t)' = \boldsymbol F(\boldsymbol u(t))$, starting with $\boldsymbol u(t_0)$ defined by the initial condition and integrating up to $t=t_1$. 

## Handling Boundary Conditions 

The above derivation assumes that $b_0, b_1 \neq 0$ and that a Robin boundary condition assumes. The `BoundaryConditions` struct can take three types of boundary conditions:

- `Dirichlet`
- `Neumann`
- `Robin`

A Dirichlet boundary condition is given by $u(a, t) = f\left(u(a, t)\right)$, and similarly for $x=b$, and so we cannot make a definition for the $a_j$ or $b_j$ coefficients in this case, with $j \in \{0, 1\}$. We instead use the callback interface from DifferentialEquations.jl for this case.

Neumann boundary conditions take the form $\partial u(a, t)/\partial x = f\left(u(a, t)\right)$, and similarly for $x=b$. This boundary condition can be written

```math
-f\left(u(a, t)\right) + \dfrac{\partial u(a, t)}{\partial x} = 0,
```

which is a Robin boundary condition with $a_0 = -f$ and $b_0 = 1$. 

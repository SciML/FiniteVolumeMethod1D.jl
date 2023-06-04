var documenterSearchIndex = {"docs":
[{"location":"examples/","page":"Examples","title":"Examples","text":"CurrentModule = FiniteVolumeMethod1D","category":"page"},{"location":"examples/#Examples","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"This section gives some examples for how this package can be used. In most of the examples they follow, there are exact solutions, but we do not discuss them here. You can see the scripts in the tests if you are interested.","category":"page"},{"location":"examples/#Example-I:-Heat-equation","page":"Examples","title":"Example I: Heat equation","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"We start with a simple example:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"beginalign*\ndfracpartial upartial t = dfracpartial^2 upartial x^2 quad 0  x  1t0 8pt \ndfracpartial u(0 t)partial x = 0 quad t0 8pt\ndfracpartial u(1 t)partial x = 0 quad t0 8pt\nu(x 0)  = x quad 0  x  1\nendalign*","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"The first step is to define the geometry and the boundary conditions. The geometry for these types of problems simply requires a set of mesh_points, which can be regularly or irregular spaced. Here we use","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"mesh_points = LinRange(0, 1, 500)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"We could then use geo = FVMGeometry(mesh_points), but we will use the simpler constructor for the FVMProblem later. Note that the constructor will take a=0 and b=1 from mesh_points[begin] and mesh_points[end]. ","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"The boundary conditions are defined using the Neumann type. Since the boundary condition is constant in this case, we can use the simpler Neumann(::Number) constructor (see ?Neumann for other constructors).","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using FiniteVolumeMethod1D\nlhs = Neumann(0.0)\nrhs = Neumann(0.0)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"As before, we could then use BoundaryConditions(lhs, rhs), but we will use the simpler constructor for the FVMProblem. Now, let us define the diffusion function, initial condition, and final time.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"diffusion_function = (u, x, t, p) -> one(u)\ninitial_condition = collect(mesh_points)\nfinal_time = 0.1","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"We can now construct the FVMProblem.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"prob = FVMProblem(mesh_points, lhs, rhs;\n    diffusion_function,\n    initial_condition,\n    final_time)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"This prob can be solved the same way as you would e.g. with DifferentialEquations.jl with solve. Using solve returns the same struct as DifferentialEquations.jl returns:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using OrdinaryDiffEq\nusing LinearSolve\nsol = solve(prob, TRBDF2(linsolve=KLUFactorization()))","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"This can be easily plotted, e.g.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using CairoMakie\nlet t_range = LinRange(0.0, final_time, 250)\n    fig = Figure(fontsize=33)\n    ax = Axis(fig[1, 1], xlabel=L\"x\", ylabel=L\"t\")\n    sol_u = [sol(t) for t in t_range]\n    contourf!(ax, mesh_points, t_range, reduce(hcat, sol_u), colormap=:viridis)\n    fig\nend","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"<figure>\n    <img src='../figures/heat_contour.png', alt'Contour of the heat equation solution'><br>\n</figure>","category":"page"},{"location":"examples/#Example-II:-Reaction-diffusion-equation","page":"Examples","title":"Example II: Reaction-diffusion equation","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"The next example we consider is a reaction-diffusion equation:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"beginalign*\ndfracpartial upartial t = dfracpartial^2 upartial x^2 + pi^2 exp-24pi^2tsin(5pi x) quad 0  x  1 t  0 8pt\nu(0 t) = 0 8pt\nu(1 t) = 0 8pt\nu(x 0) = 3sin(4pi x)\nendalign*","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"The boundary conditions can be constructed using Dirichlet. So, the geometry and boundary conditions can be given by:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using FiniteVolumeMethod1D\nmesh_points = LinRange(0, 1, 500)\nlhs = Dirichlet(0.0)\nrhs = Dirichlet(0.0)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Now, let us define the diffusion and reaction functions. We show here the use of parameters.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"diffusion_function = (u, x, t, p) -> one(u)\nreaction_function = (u, x, t, p) -> p[1] * exp(-p[2] * t) * sin(p[3] * x)\nreaction_parameters = (π^2, 24π^2, 5π)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Now the problem can be constructed.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"ic_ff = x -> 3sin(4π*x)\ninitial_condition = ic_ff.(mesh_points)\nfinal_time = 0.05\nprob = FVMProblem(mesh_points, lhs, rhs;\n    diffusion_function,\n    reaction_function,\n    reaction_parameters,\n    initial_condition,\n    final_time)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Now we solve and plot the solution.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using OrdinaryDiffEq \nusing LinearSolve \nsol = solve(prob, TRBDF2(linsolve=KLUFactorization()), saveat = 0.001)\n\nusing CairoMakie\nlet t_range = LinRange(0.0, final_time, 250)\n    fig = Figure(fontsize=33)\n    ax = Axis3(fig[1, 1], xlabel=L\"x\", ylabel=L\"t\", zlabel=L\"z\", azimuth = 0.8)\n    sol_u = [sol(t) for t in t_range]\n    surface!(ax, mesh_points, t_range, reduce(hcat, sol_u), colormap=:viridis)\n    fig\nend","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"<figure>\n    <img src='../figures/dirichlet_source_surface.png', alt'Surface plot of the reaction-diffusion problem'><br>\n</figure>","category":"page"},{"location":"examples/#Example-III:-Diffusion-problem-with-Robin-boundary-conditions","page":"Examples","title":"Example III: Diffusion problem with Robin boundary conditions","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"This next problem we consider is a diffusion problem with Robin boundary conditions:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"beginalign*\ndfracpartial upartial t  = frac125dfracpartial^2 upartial x^2 quad 0  x  3 t0 8pt\nu(0 t)  = 0 8pt\ndfrac12u(3 t) + dfracpartial u(3 t)partial x = 08pt \nu(x 0)  = 100left(1-dfracx3right)\nendalign*","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Robin boundary conditions are supported by rewriting them in Neumann form, so that partial u(3 t)partial x = -u(3 t)2 in the above. Thus:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using FiniteVolumeMethod1D\nmesh_points = LinRange(0, 3, 2500)\nlhs = Dirichlet(0.0)\nrhs = Neumann((u, t, p) -> u * p, -1 / 2)\ndiffusion_function = (u, x, t, p) -> p^2\ndiffusion_parameters = 1 / 5\nic = x -> 100(1 - x / 3)\ninitial_condition = ic.(mesh_points)\nfinal_time = 3.0\nprob = FVMProblem(mesh_points, lhs, rhs;\n    diffusion_function,\n    diffusion_parameters,\n    final_time=final_time,\n    initial_condition\n)\n\nusing OrdinaryDiffEq \nusing LinearSolve\nsol = solve(prob, TRBDF2(linsolve=KLUFactorization()), saveat = 0.01)\n\nusing CairoMakie\nlet t_range = LinRange(0.0, final_time, 250)\n    fig = Figure(fontsize=33)\n    ax = Axis3(fig[1, 1], xlabel=L\"x\", ylabel=L\"t\", zlabel=L\"z\", azimuth = 0.8)\n    sol_u = [sol(t) for t in t_range]\n    surface!(ax, mesh_points, t_range, reduce(hcat, sol_u), colormap=:viridis)\n    fig\nend","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"<figure>\n    <img src='../figures/robin_diffusion_surface.png', alt'Surface plot of the Robin diffusion problem'><br>\n</figure>","category":"page"},{"location":"examples/#Example-IV:-Porous-Fisher-equation-with-degenerate-diffusion","page":"Examples","title":"Example IV: Porous-Fisher equation with degenerate diffusion","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"We now consider the following problem:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"beginalign*\ndfracpartial upartial t = dfracpartialpartial xleft(D(u)dfracpartial upartial xright) + R(u) quad -6pi  x  6pi8pt \ndfracpartial u(-2pi t)partial x  = 0 8pt\ndfracpartial u(2pi t)partial x  = 0 8pt \nu(x 0)  = begincases 1  x  0  12  x geq 0 endcases\nendalign*","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"where D(u) = 1(10u) + 50u^2 + 3u^3 and R(u) = beta K u(1 - uK). We take beta = 10^-3 and K = 2. The problem is solved as follows:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using FiniteVolumeMethod1D \nmesh_points = LinRange(-2π, 2π, 500)\nlhs = Neumann(0.0)\nrhs = Neumann(0.0)\ndiffusion_function = (u, x, t, p) -> inv(p[1] * u) + p[2] * inv(u^2) + p[3] * inv(u^3)\ndiffusion_parameters = [1.0, 50.0, 3.0]\nreaction_function = (u, x, t, p) -> p[1] * p[2] * u * (1 - u / p[2])\nreaction_parameters = [1e-3, 2.0]\nic_f(x) = x < 0 ? 1.0 : 1 / 2\ninitial_condition = ic_f.(mesh_points)\nfinal_time = 1.0\nprob = FVMProblem(\n    mesh_points,\n    lhs,\n    rhs;\n    diffusion_function,\n    reaction_function,\n    diffusion_parameters,\n    reaction_parameters,\n    initial_condition,\n    final_time\n)\n\nusing OrdinaryDiffEq \nsol = solve(fvm_prob, TRBDF2())\n\nusing CairoMakie \nlet t_range = LinRange(0.0, final_time, 250)\n    fig = Figure(fontsize=33)\n    ax = Axis3(fig[1, 1], xlabel=L\"x\", ylabel=L\"t\", zlabel=L\"z\", azimuth = 0.8)\n    sol_u = [sol(t) for t in t_range]\n    surface!(ax, mesh_points, t_range, reduce(hcat, sol_u), colormap=:viridis)\n    fig\nend","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"<figure>\n    <img src='../figures/fisher_surface.png', alt'Surface plot of the Porous-Fisher problem'><br>\n</figure>","category":"page"},{"location":"math/","page":"Mathematical Details","title":"Mathematical Details","text":"CurrentModule = FiniteVolumeMethod1D","category":"page"},{"location":"math/#Mathematical-Details","page":"Mathematical Details","title":"Mathematical Details","text":"","category":"section"},{"location":"math/","page":"Mathematical Details","title":"Mathematical Details","text":"In this section, we provide some of the mathematical details for discretising the PDEs we consider. Recall that the problems we consider are","category":"page"},{"location":"math/","page":"Mathematical Details","title":"Mathematical Details","text":"beginalign*\ndfracpartial upartial t = dfracpartialpartial xleft(D(u)dfracpartial upartial xright) + R(u)  a leq x leq b t_0  t leq t_1 9pt\ndfracpartial u(a t)partial x = a_0left(u(a t) tright)  t_0  t leq t_1 9pt\ndfracpartial u(b t)partial x = a_1left(u(b t) tright)  t_0  t leq t_1 9pt\nu(x 0) = f(x)  a leq x leq b\nendalign*","category":"page"},{"location":"math/","page":"Mathematical Details","title":"Mathematical Details","text":"This is for the Neumann boundary condition form at both ends. Dirichlet boundary conditions are handled via callbacks, as discussed in the examples; we assume that b_0 b_1 neq 0 in what follows. (We also support functions with arguments x and t, e.g. D(u x t), but for simplicity we omit the x and t arguments.)","category":"page"},{"location":"math/#Interior-Discretisation","page":"Mathematical Details","title":"Interior Discretisation","text":"","category":"section"},{"location":"math/","page":"Mathematical Details","title":"Mathematical Details","text":"Let us start by focusing on the discretisation of the PDE itself. We start by defining some grid x_1 ldots x_n for the mesh points, assuming a = x_1  x_2  cdots  x_n = b. The control volumes are defined by intervals w_i e_i, where","category":"page"},{"location":"math/","page":"Mathematical Details","title":"Mathematical Details","text":"beginalign*\nw_i = begincases x_1  i=1  dfrac12left(x_i-1 + x_iright)  i=2ldotsn endcases \ne_i = begincases frac12left(x_i + x_i+1right)  i=1ldotsn-1  x_n  i=n endcases \nendalign*","category":"page"},{"location":"math/","page":"Mathematical Details","title":"Mathematical Details","text":"The volumes of these control volumes are defined by V_i = e_i - w_i, i=1ldotsn. Now, integrate the PDE over a control volume w_i e_i:","category":"page"},{"location":"math/","page":"Mathematical Details","title":"Mathematical Details","text":"beginalign*\nint_w_i^e_i dfracpartial upartial tmathrmdx = int_w_i^e_i dfracpartialpartial xleft(D(u)dfracpartial upartial xright)mathrmdx + int_w_i^e_i R(u)mathrmdx \ndfracmathrm dmathrm dtint_w_i^e_i umathrmdx = Dleft(u(e_i t)right)dfracpartial u(e_i t)partial x - Dleft(u(w_i t)right)dfracpartial u(w_i t)partial x + int_w_i^e_i R(u)mathrmdx \ndfracmathrm dbar u_imathrm dt = frac1V_ileftDleft(u(e_i t)right)dfracpartial u(e_i t)partial x - Dleft(u(w_i t)right)dfracpartial u(w_i t)partial xright + bar R_i\nendalign*","category":"page"},{"location":"math/","page":"Mathematical Details","title":"Mathematical Details","text":"where bar u_i = (1V_i)int_w_i^e_i umathrmdx and bar R_i = (1V_i)int_w_i^e_i R(u)mathrmdx. Letting u_i = u(x_i t) and R_i = R(u_i), we make the following approximations:","category":"page"},{"location":"math/","page":"Mathematical Details","title":"Mathematical Details","text":"beginalign*\nbeginarrayrcll\nbar u_i approx u_i  i=1ldotsn 6pt\nbar R_i approx R_i  i=1ldotsn 6pt\nDleft(u(e_i t)right) approx dfrac12left(D_i + D_i+1right)quad  i=1ldotsn-1 6pt\nDleft(u(w_i t)right) approx dfrac12left(D_i-1 + D_iright)quad  i=2ldotsn 6pt\ndfracpartial u(e_i t)partial x approx dfracu_i+1 - u_ih_i  i=1ldotsn-1 6pt\ndfracpartial u(w_i t)partial x approx dfracu_i - u_i-1h_i-1  i=2ldotsn \nendarray\nendalign*","category":"page"},{"location":"math/","page":"Mathematical Details","title":"Mathematical Details","text":"where h_i = x_i+1 - x_i, i=1ldotsn-1. With these approximations, we find:","category":"page"},{"location":"math/","page":"Mathematical Details","title":"Mathematical Details","text":"beginalign*\nfracmathrm du_imathrm dt = frac1V_ileftleft(dfracD_i+D_i+12right)left(dfracu_i+1 - u_ih_iright) - left(dfracD_i-1 + D_i2right)left(dfracu_i - u_i-1h_i-1right)right + R_i\nendalign*","category":"page"},{"location":"math/","page":"Mathematical Details","title":"Mathematical Details","text":"for i=2ldotsn-1.","category":"page"},{"location":"math/#Boundary-Discretisation","page":"Mathematical Details","title":"Boundary Discretisation","text":"","category":"section"},{"location":"math/","page":"Mathematical Details","title":"Mathematical Details","text":"We still need to handle the equations at i=1 and i=n. Using partial u(a t) = a_0(u(a t) t), we obtain Thus,","category":"page"},{"location":"math/","page":"Mathematical Details","title":"Mathematical Details","text":"dfracmathrm du_1mathrm dt = frac1V_1leftleft(dfracD_1 + D_22right)left(dfracu_2 - u_1h_1right) - a_0(u_1 t)D(u_1)right + R_1","category":"page"},{"location":"math/","page":"Mathematical Details","title":"Mathematical Details","text":"Similarly, ","category":"page"},{"location":"math/","page":"Mathematical Details","title":"Mathematical Details","text":"dfracmathrm du_nmathrm dt = frac1V_nlefta_1(u_1)D(u_n) - left(dfracD_n-1 + D_n2right)left(dfracu_n - u_n-1h_n-1right)right + R_n","category":"page"},{"location":"math/#The-Complete-Discretisation","page":"Mathematical Details","title":"The Complete Discretisation","text":"","category":"section"},{"location":"math/","page":"Mathematical Details","title":"Mathematical Details","text":"Putting all the results together, the complete system of ODEs is","category":"page"},{"location":"math/","page":"Mathematical Details","title":"Mathematical Details","text":"beginalign*\nfracmathrm du_imathrm dt = frac1V_ileftleft(dfracD_i+D_i+12right)left(dfracu_i+1 - u_ih_iright) - left(dfracD_i-1 + D_i2right)left(dfracu_i - u_i-1h_i-1right)right + R_i i=2ldotsn-1 8pt\ndfracmathrm du_1mathrm dt = frac1V_1leftleft(dfracD_1 + D_22right)left(dfracu_2 - u_1h_1right) - a_0(u_1 t)D(u_1)right + R_18pt\ndfracmathrm du_nmathrm dt=frac1V_nlefta_1(u_1)D(u_n) - left(dfracD_n-1 + D_n2right)left(dfracu_n - u_n-1h_n-1right)right + R_n\nendalign*","category":"page"},{"location":"math/","page":"Mathematical Details","title":"Mathematical Details","text":"This system can then be easily solved using methods from DifferentialEquation.jl, treating the system in the form boldsymbol u(t) = boldsymbol F(boldsymbol u(t)), starting with boldsymbol u(t_0) defined by the initial condition and integrating up to t=t_1. ","category":"page"},{"location":"math/#Handling-Boundary-Conditions","page":"Mathematical Details","title":"Handling Boundary Conditions","text":"","category":"section"},{"location":"math/","page":"Mathematical Details","title":"Mathematical Details","text":"The above derivation assumes that b_0 b_1 neq 0 and that a Neumann boundary condition is used. The BoundaryConditions struct can take two types of boundary conditions:","category":"page"},{"location":"math/","page":"Mathematical Details","title":"Mathematical Details","text":"Dirichlet\nNeumann","category":"page"},{"location":"math/","page":"Mathematical Details","title":"Mathematical Details","text":"A Dirichlet boundary condition is given by u(a t) = gleft(u(a t) tright), and similarly for x=b, and so we cannot make a definition for the a_j or b_j coefficients in this case, with j in 0 1. We instead use the callback interface from DifferentialEquations.jl for this case.","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = FiniteVolumeMethod1D","category":"page"},{"location":"#FiniteVolumeMethod1D","page":"Home","title":"FiniteVolumeMethod1D","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for FiniteVolumeMethod1D. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"This is a package for solving equations of the form","category":"page"},{"location":"","page":"Home","title":"Home","text":"dfracpartial u(x t)partial x = dfracpartialpartial xleft(Dleft(u x tright)dfracpartial u(x t)partial xright) + R(u x t)","category":"page"},{"location":"","page":"Home","title":"Home","text":"using the finite volume method over intervals a leq x leq b and t_0 leq t leq t_1, with support for the following types of boundary conditions (shown at x = a, but you can mix boundary condition types, e.g. Neumann at x=a and Robin at x=b):","category":"page"},{"location":"","page":"Home","title":"Home","text":"beginalign*\nbeginarrayrrcl\ntextRobin  a_0left(u(a t) tright) + b_0left(u(a t) tright)dfracpartial u(a t)partial t  =  0 9pt\ntextNeumann  dfracpartial u(a t)partial t  =  a_0left(u(a t) tright) 9pt\ntextDirichlet  u(a t)  =  a_0left(u(a t) tright)\nendarray\nendalign*","category":"page"},{"location":"","page":"Home","title":"Home","text":"where the Dirichlet condition has u(a t) mapping from a_0(u(a t) t) (i.e., it is not an implicit equation for u(a t)).","category":"page"},{"location":"","page":"Home","title":"Home","text":"The package is not registered, so to install it you must do:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> ] add https://github.com/DanielVandH/FiniteVolumeMethod1D.jl\njulia> using FiniteVolumeMethod1D","category":"page"},{"location":"","page":"Home","title":"Home","text":"More information is given in the sidebar, and the docstrings are below.","category":"page"},{"location":"","page":"Home","title":"Home","text":"If you want a more complete two-dimensional version, please see my other (registered) package FiniteVolumeMethod.jl.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [FiniteVolumeMethod1D]","category":"page"},{"location":"#FiniteVolumeMethod1D.BoundaryConditions","page":"Home","title":"FiniteVolumeMethod1D.BoundaryConditions","text":"BoundaryConditions{L<:AbstractBoundaryCondition,R<:AbstractBoundaryCondition}\n\nThe boundary conditions for the FVMProblem. \n\nFields\n\nlhs::L: The left-hand side boundary condition.\nrhs::R: The right-hand side boundary condition.\n\nSee also Dirichlet and Neumann for the types of  boundary conditions you can construct.\n\n\n\n\n\n","category":"type"},{"location":"#FiniteVolumeMethod1D.Dirichlet","page":"Home","title":"FiniteVolumeMethod1D.Dirichlet","text":"Dirichlet{F,P} <: AbstractBoundaryCondition{F,P}\n\nA Dirichlet boundary condition with fields f and p (default p = nothing), with f being a function of the form f(u, p) and p being the parameters for f. \n\nA Dirichlet boundary condition takes the form\n\nu(a t)  f(u(a t) t p)\n\nwhere a is one of the endpoints. \n\nConstructors\n\nDirichlet(f::Function, p = nothing) -> Dirichlet(f, p)\nDirichlet(; f, p = nothing)         -> Dirichlet(f, p)\nDirichlet(v::Number)                -> Dirichlet((u, t, p) -> oftype(u, v), nothing)\n\n\n\n\n\n","category":"type"},{"location":"#FiniteVolumeMethod1D.FVMGeometry","page":"Home","title":"FiniteVolumeMethod1D.FVMGeometry","text":"FVMGeometry{T}\n\nDefinition of the geometry for a finite volume method problem.\n\nFields\n\nmesh_points::T: The mesh points. Must be sorted.\nspacings::T: The spacings between the mesh points. \nvolumes::T: The volumes of the cells defined by the mesh points.\n\nConstructors\n\nTo construct the geometry, you can directly call the default constructor, \n\nFVMGeometry(mesh_points, spacings, volumes)\n\nor you can call the convenience constructor,\n\nFVMGeometry(mesh_points)\n\nwhich will compute the spacings and volumes for you.\n\nSee also FVMProblem.\n\n\n\n\n\n","category":"type"},{"location":"#FiniteVolumeMethod1D.FVMProblem","page":"Home","title":"FiniteVolumeMethod1D.FVMProblem","text":"FVMProblem{T,DF,DP,RF,RP,L,R,IC,FT}\n\nDefinition of an FVMProblem.\n\nFields\n\ngeometry::FVMGeometry{T}: The geometry of the problem.\nboundary_conditions::BoundaryConditions{L, R}: The boundary conditions.\ndiffusion_function::DF: The diffusion function, of the form (u, x, t, p) -> Number, where p = diffusion_parameters.\ndiffusion_parameters::DP: The parameters for the diffusion function.\nreaction_function::RF: The reaction function, of the form (u, x, t, p) -> Number, where p = reaction_parameters.\nreaction_parameters::RP: The parameters for the reaction function.\ninitial_condition::IC: The initial condition, with initial_condition[i] corresponding to the value at geometry.mesh_points[i] and t = initial_time.\ninitial_time::FT: The initial time.\nfinal_time::FT: The final time.\n\nConstructors\n\nYou can use the default constructor, but we also provide the constructor \n\nFVMProblem(;\n    geometry, \n    boundary_conditions,\n    diffusion_function,\n    diffusion_parameters = nothing,\n    reaction_function = Returns(0.0),\n    reaction_parameters = nothing,\n    initial_condition,\n    initial_time = 0.0,\n    final_time)\n\nwhich provides some default values. Moreover, instead of providing geometry and boundary_conditions, you can use \n\nFVMProblem(mesh_points, lhs, lhs; kwargs...)\n\nwhich will construct geometry = FVMGeometry(mesh_points) and boundary_conditions = BoundaryConditions(lhs, rhs).  The kwargs... are as above, except without geometry and boundary_conditions of course.\n\nTo solve the FVMProblem, just use solve as you would in DifferentialEquations.jl. For example, \n\nsol = solve(prob, Tsit5(), saveat=0.1)\n\n\n\n\n\n","category":"type"},{"location":"#FiniteVolumeMethod1D.Neumann","page":"Home","title":"FiniteVolumeMethod1D.Neumann","text":"Neumann{F,P} <: AbstractBoundaryCondition{F,P}\n\nA Neumann boundary condition with fields f and p (default p = nothing), with f being a function of the form f(u, t, p) and p being the parameters for f.\n\nA Neumann boundary condition takes the form\n\ndfracpartial upartial x(a t) = f(u(a t) t p)\n\nwhere a is one of the endpoints. \n\nConstructors\n\nNeumann(f::Function, p = nothing) -> Neumann(f, p)\nNeumann(; f, p = nothing)         -> Neumann(f, p)\nNeumann(v::Number)                -> Neumann((u, t, p) -> oftype(u, v), nothing)\n\n\n\n\n\n","category":"type"}]
}

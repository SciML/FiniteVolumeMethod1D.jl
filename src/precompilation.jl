using PrecompileTools

@setup_workload begin
    # Minimal setup - these are lightweight and don't add much to precompilation time
    mesh_points = LinRange(0.0, 1.0, 10)

    @compile_workload begin
        # Precompile FVMGeometry creation - common entry point
        geom = FVMGeometry(mesh_points)

        # Precompile boundary condition types with common cases
        lhs_dir = Dirichlet(0.0)
        rhs_dir = Dirichlet(1.0)
        lhs_neu = Neumann(0.0)
        rhs_neu = Neumann(0.0)

        # Precompile BoundaryConditions with different combinations
        bc_dir = BoundaryConditions(lhs = lhs_dir, rhs = rhs_dir)
        bc_neu = BoundaryConditions(lhs = lhs_neu, rhs = rhs_neu)
        bc_mix = BoundaryConditions(lhs = lhs_dir, rhs = rhs_neu)

        # Precompile FVMProblem creation
        diffusion_function = (u, x, t, p) -> one(u)
        initial_condition = collect(mesh_points)

        prob = FVMProblem(
            mesh_points, lhs_dir, rhs_dir;
            diffusion_function = diffusion_function,
            initial_condition = initial_condition,
            final_time = 1.0
        )

        # Precompile ODEProblem creation - this is a common operation
        ode_prob = SciMLBase.ODEProblem(prob)

        # Precompile jacobian_sparsity
        jac_sp = jacobian_sparsity(prob)
    end
end

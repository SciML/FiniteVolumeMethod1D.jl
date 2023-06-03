function SciMLBase.ODEProblem(prob::FVMProblem;
    specialization::Type{S}=SciMLBase.AutoSpecialize,
    jac_prototype=jacobian_sparsity(prob),
    kwargs...) where {S}
    initial_time = prob.initial_time
    final_time = prob.final_time
    time_span = (initial_time, final_time)
    initial_condition = prob.initial_condition
    kwarg_dict = Dict(kwargs)
    lhs_cb = lhs_dirichlet(prob, :saveat ∈ keys(kwarg_dict))
    rhs_cb = rhs_dirichlet(prob, :saveat ∈ keys(kwarg_dict))
    cb = CallbackSet(lhs_cb, rhs_cb)
    f = ODEFunction{true,S}(pde_odes!; jac_prototype=jac_prototype)
    ode_problem = ODEProblem{true,S}(f, initial_condition, time_span, prob; callback=cb, kwargs...)
    return ode_problem
end

@inline function lhs_dirichlet(prob, has_saveat=true)
    cb = let is_dir = is_dirichlet(prob.boundary_conditions.lhs), has_sv = has_saveat
        condition = (u, t, integrator) -> is_dir
        DiscreteCallback(
            condition,
            update_lhs!;
            save_positions=(!has_sv, !has_sv)
        )
    end
    return cb
end
@inline function rhs_dirichlet(prob, has_saveat=true)
    cb = let is_dir = is_dirichlet(prob.boundary_conditions.rhs), has_sv = has_saveat
        condition = (u, t, integrator) -> is_dir
        DiscreteCallback(
            condition,
            update_rhs!;
            save_positions=(!has_sv, !has_sv)
        )
    end
    return cb
end
function update_lhs!(integrator)
    u = integrator.u
    t = integrator.t
    prob = integrator.p
    boundary_conditions = prob.boundary_conditions
    lhs = boundary_conditions.lhs
    val = lhs(u[begin], t)
    u[begin] = val
    return nothing
end
function update_rhs!(integrator)
    u = integrator.u
    t = integrator.t
    prob = integrator.p
    boundary_conditions = prob.boundary_conditions
    rhs = boundary_conditions.rhs
    val = rhs(u[begin], t)
    u[end] = val
    return nothing
end

CommonSolve.init(prob::FVMProblem, alg; kwargs...) = CommonSolve.init(ODEProblem(prob; kwargs...), alg; kwargs...)

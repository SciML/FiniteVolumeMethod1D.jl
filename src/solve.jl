function SciMLBase.ODEProblem(prob::FVMProblem;
    specialization::Type{S}=SciMLBase.AutoSpecialize,
    jac_prototype=jacobian_sparsity(prob),
    kwargs...) where {S}
    initial_time = prob.initial_time
    final_time = prob.final_time
    time_span = (initial_time, final_time)
    initial_condition = prob.initial_condition
    f = ODEFunction{true,S}(pde_odes!; jac_prototype=jac_prototype)
    ode_problem = ODEProblem{true,S}(f, initial_condition, time_span, prob; kwargs...)
    return ode_problem
end

CommonSolve.init(prob::FVMProblem, alg; kwargs...) = CommonSolve.init(ODEProblem(prob; kwargs...), alg; kwargs...)

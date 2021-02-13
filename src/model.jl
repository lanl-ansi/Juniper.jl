include("debug.jl")
include("fpump.jl")

function create_root_model!(optimizer::MOI.AbstractOptimizer, jp::JuniperProblem; fix_start=false)
    ps = jp.options.log_levels

    jp.model = MOI.instantiate(jp.nl_solver, with_bridge_type=Float64)
    # all continuous we solve relaxation first
    jp.x = Vector{MOI.VariableIndex}(undef, jp.num_var)
    jp.cx = Vector{MOI.ConstraintIndex{MOI.SingleVariable, MOI.Interval{Float64}}}(undef, jp.num_var)
    for i in 1:jp.num_var
        jp.x[i], jp.cx[i] = MOI.add_constrained_variable(jp.model, MOI.Interval(jp.l_var[i], jp.u_var[i]))
    end
    # FIXME we should use `copy_to` here to actually map indices and support
    # cases where indices of the optimizer are not as expected.
    for i in eachindex(jp.x)
        @assert jp.x[i] == MOI.VariableIndex(i)
    end
    MOI.set(jp.model, MOI.VariablePrimalStart(), jp.x, jp.primal_start)

    MOI.set(jp.model, MOI.ObjectiveSense(), optimizer.sense)
    MOI.set(jp.model, MOI.NLPBlock(), optimizer.nlp_data)

    # Nonlinear objectives *override* any objective set by using the `ObjectiveFunction` attribute.
    # So we first check that `optimizer.nlp_data.has_objective` is `false`.
    if !optimizer.nlp_data.has_objective && optimizer.objective !== nothing
        MOI.set(jp.model, MOI.ObjectiveFunction{typeof(optimizer.objective)}(), optimizer.objective)
    end

    llc = optimizer.linear_le_constraints
    lgc = optimizer.linear_ge_constraints
    lec = optimizer.linear_eq_constraints
    qlc = optimizer.quadratic_le_constraints
    qgc = optimizer.quadratic_ge_constraints
    qec = optimizer.quadratic_eq_constraints
    for constr_type in [llc, lgc, lec, qlc, qgc, qec]
        for constr in constr_type
            MOI.add_constraint(jp.model, constr[1], constr[2])
        end
    end
end

function fix_primal_start!(jp::JuniperProblem)
    lb = jp.l_var
    ub = jp.u_var
    x = jp.x
    for i=1:jp.num_var
        if jp.var_type[i] != :Cont && lb[i] <= jp.primal_start[i] <= ub[i]
            MOI.set(jp.model, MOI.ConstraintSet(), jp.cx[i], MOI.Interval(jp.primal_start[i], jp.primal_start[i]))
        else
            MOI.set(jp.model, MOI.VariablePrimalStart(), jp.x[i], jp.primal_start[i])
        end
    end
end

function unfix_primal_start!(jp::JuniperProblem)
    lb = jp.l_var
    ub = jp.u_var
    x = jp.x
    for i=1:jp.num_var
        MOI.set(jp.model, MOI.ConstraintSet(), jp.cx[i], MOI.Interval(lb[i], ub[i]))
        MOI.set(jp.model, MOI.VariablePrimalStart(), jp.x[i], jp.primal_start[i])
    end
end

function solve_root_incumbent_model(jp::JuniperProblem)
    MOI.optimize!(jp.model)
    status = MOI.get(jp.model, MOI.TerminationStatus())
    incumbent = nothing
    ps = jp.options.log_levels
    if state_is_optimal(status; allow_almost=jp.options.allow_almost_solved)
        # set incumbent
        objval = MOI.get(jp.model, MOI.ObjectiveValue())
        solution = MOI.get(jp.model, MOI.VariablePrimal(), jp.x)
        if are_type_correct(solution, jp.var_type, jp.disc2var_idx, jp.options.atol)
            if only_almost_solved(status) && jp.options.allow_almost_solved_integral
                @warn "Start value incumbent only almost locally solved. Disable with `allow_almost_solved_integral=false`"
            end

            incumbent = Incumbent(objval, solution, only_almost_solved(status))
            (:All in ps || :Info in ps) && println("Incumbent using start values: $objval")
        end
    else
        (:All in ps || :Info in ps) && println("Start values are not feasible.")
    end
    return incumbent
end

function solve_root_model!(jp::JuniperProblem)
    MOI.optimize!(jp.model)
    status = MOI.get(jp.model, MOI.TerminationStatus())
    jp.relaxation_status = status
    restarts = 0
    max_restarts = jp.options.num_resolve_root_relaxation
    jp.options.debug && debug_init(jp.debugDict)
    while !state_is_optimal(jp.relaxation_status; allow_almost=jp.options.allow_almost_solved) &&
        restarts < max_restarts && time()-jp.start_time < jp.options.time_limit

        # TODO freemodel for Knitro
        restart_values = generate_random_restart(jp)
        jp.options.debug && debug_restart_values(jp.debugDict,restart_values)
        MOI.set(jp.model, MOI.VariablePrimalStart(), jp.x, restart_values)
        MOI.optimize!(jp.model)
        jp.relaxation_status = MOI.get(jp.model, MOI.TerminationStatus())
        restarts += 1
    end
    if only_almost_solved(jp.relaxation_status)
        @warn "The relaxation is only almost solved."
    end
    return restarts
end

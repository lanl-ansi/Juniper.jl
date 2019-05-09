include("debug.jl")
include("fpump.jl")

function create_root_model!(optimizer::MOI.AbstractOptimizer, jp::JuniperProblem)
    ps = jp.options.log_levels

    jp.model = Model()
    lb = jp.l_var
    ub = jp.u_var
    # all continuous we solve relaxation first
    @variable(jp.model, lb[i] <= x[i=1:jp.num_var] <= ub[i])

    for i=1:jp.num_var
        JuMP.set_start_value(x[i], jp.primal_start[i])
    end

    if jp.options.registered_functions != nothing
        for reg_f in jp.options.registered_functions
            if reg_f.gradf == nothing
                JuMP.register(jp.model, reg_f.s, reg_f.dimension, reg_f.f; autodiff=reg_f.autodiff)
            elseif reg_f.grad2f == nothing
                JuMP.register(jp.model, reg_f.s, reg_f.dimension, reg_f.f, reg_f.gradf)
            else
                JuMP.register(jp.model, reg_f.s, reg_f.dimension, reg_f.f, reg_f.gradf, reg_f.grad2f)
            end
        end
    end

    # TODO check whether it is supported
    if optimizer.nlp_data.has_objective
        obj_expr = MOI.objective_expr(optimizer.nlp_data.evaluator)
        expr_dereferencing!(obj_expr, jp.model)
        JuMP.set_NL_objective(jp.model, optimizer.sense, obj_expr)
    elseif optimizer.objective != nothing
        MOI.set(jp.model, MOI.ObjectiveFunction{typeof(optimizer.objective)}(), optimizer.objective)
        MOI.set(jp.model, MOI.ObjectiveSense(), optimizer.sense)
    end

    backend = JuMP.backend(jp.model);
    llc = optimizer.linear_le_constraints
    lgc = optimizer.linear_ge_constraints
    lec = optimizer.linear_eq_constraints
    qlc = optimizer.quadratic_le_constraints
    qgc = optimizer.quadratic_ge_constraints
    qec = optimizer.quadratic_eq_constraints
    for constr_type in [llc, lgc, lec, qlc, qgc, qec]
        for constr in constr_type
            MOI.add_constraint(backend, constr[1], constr[2])
        end
    end
    for i in 1:jp.num_nl_constr
        constr_expr = MOI.constraint_expr(optimizer.nlp_data.evaluator, i)
        expr_dereferencing!(constr_expr, jp.model)
        JuMP.add_NL_constraint(jp.model, constr_expr)
    end

    (:All in ps || :Info in ps) && print_info(jp)
    
    jp.x = x
end

function solve_root_model!(jp::JuniperProblem)
    status, backend = optimize_get_status_backend(jp.model; solver=jp.nl_solver)
    jp.status = status
    restarts = 0
    max_restarts = jp.options.num_resolve_root_relaxation
    jp.options.debug && debug_init(jp.debugDict)
    while !state_is_optimal(jp.status; allow_almost=true) &&
        restarts < max_restarts && time()-jp.start_time < jp.options.time_limit

        # TODO freemodel for Knitro
        restart_values = generate_random_restart(jp)
        jp.options.debug && debug_restart_values(jp.debugDict,restart_values)
        for i=1:jp.num_var
            JuMP.set_start_value(jp.x[i], restart_values[i])
        end
        JuMP.optimize!(jp.model)
        jp.status = MOI.get(backend, MOI.TerminationStatus()) 
        restarts += 1
    end
    if jp.status == MOI.ALMOST_LOCALLY_SOLVED
        @warn "The relaxation is only almost locally solved."
    end
    return restarts
end
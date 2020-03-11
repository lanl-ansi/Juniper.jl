include("debug.jl")
include("fpump.jl")

function create_root_model!(optimizer::MOI.AbstractOptimizer, jp::JuniperProblem; fix_start=false)
    ps = jp.options.log_levels

    jp.model = Model()
    lb = jp.l_var
    ub = jp.u_var
    # all continuous we solve relaxation first
    @variable(jp.model, lb[i] <= x[i=1:jp.num_var] <= ub[i])

    for i=1:jp.num_var
        JuMP.set_start_value(x[i], jp.primal_start[i])
    end
    
    if jp.options.registered_functions !== nothing
        for reg_f in jp.options.registered_functions
            if reg_f.gradf === nothing
                JuMP.register(jp.model, reg_f.s, reg_f.dimension, reg_f.f; autodiff=reg_f.autodiff)
            elseif reg_f.grad2f === nothing
                JuMP.register(jp.model, reg_f.s, reg_f.dimension, reg_f.f, reg_f.gradf)
            else
                JuMP.register(jp.model, reg_f.s, reg_f.dimension, reg_f.f, reg_f.gradf, reg_f.grad2f)
            end
        end
    end

    # TODO check whether it is supported
    if optimizer.nlp_data.has_objective
        obj_expr = MOI.objective_expr(optimizer.nlp_data.evaluator)
        obj_expr = expr_dereferencing(obj_expr, jp.model)
        try
            JuMP.set_NL_objective(jp.model, optimizer.sense, obj_expr)
        catch 
            error("Have you registered a function? Then please register the function also for Juniper see: https://lanl-ansi.github.io/Juniper.jl/stable/options/#registered_functions::Union{Nothing,Vector{RegisteredFunction}}-[nothing]-1")
        end
    elseif optimizer.objective !== nothing
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
        constr_expr = expr_dereferencing(constr_expr, jp.model)
        try
            JuMP.add_NL_constraint(jp.model, constr_expr)
        catch 
            error("Have you registered a function? Then please register the function also for Juniper see: https://lanl-ansi.github.io/Juniper.jl/stable/options/#registered_functions::Union{Nothing,Vector{RegisteredFunction}}-[nothing]-1")
        end
    end
    
    jp.x = x
end

function fix_primal_start!(jp::JuniperProblem)
    lb = jp.l_var
    ub = jp.u_var
    x = jp.x
    for i=1:jp.num_var
        if jp.var_type[i] != :Cont && lb[i] <= jp.primal_start[i] <= ub[i]
            JuMP.set_lower_bound(jp.x[i], jp.primal_start[i])
            JuMP.set_upper_bound(jp.x[i], jp.primal_start[i])
        else
            JuMP.set_start_value(x[i], jp.primal_start[i])
        end
    end
end

function unfix_primal_start!(jp::JuniperProblem)
    lb = jp.l_var
    ub = jp.u_var
    x = jp.x
    for i=1:jp.num_var
        JuMP.set_lower_bound(jp.x[i], lb[i])
        JuMP.set_upper_bound(jp.x[i], ub[i])
        JuMP.set_start_value(x[i], jp.primal_start[i])
    end
end

function solve_root_incumbent_model(jp::JuniperProblem)
    status, backend = optimize_get_status_backend(jp.model; solver=jp.nl_solver)
    incumbent = nothing
    ps = jp.options.log_levels
    if state_is_optimal(status; allow_almost=jp.options.allow_almost_solved)
        # set incumbent
        objval = MOI.get(backend, MOI.ObjectiveValue())
        solution = JuMP.value.(jp.x)
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
    status, backend = optimize_get_status_backend(jp.model; solver=jp.nl_solver)
    jp.relaxation_status = status
    restarts = 0
    max_restarts = jp.options.num_resolve_root_relaxation
    jp.options.debug && debug_init(jp.debugDict)
    while !state_is_optimal(jp.relaxation_status; allow_almost=jp.options.allow_almost_solved) &&
        restarts < max_restarts && time()-jp.start_time < jp.options.time_limit

        # TODO freemodel for Knitro
        restart_values = generate_random_restart(jp)
        jp.options.debug && debug_restart_values(jp.debugDict,restart_values)
        for i=1:jp.num_var
            JuMP.set_start_value(jp.x[i], restart_values[i])
        end
        JuMP.optimize!(jp.model)
        jp.relaxation_status = MOI.get(backend, MOI.TerminationStatus()) 
        restarts += 1
    end
    if only_almost_solved(jp.relaxation_status)
        @warn "The relaxation is only almost solved."
    end
    return restarts
end
include("debug.jl")
include("fpump.jl")

function create_root_model!(
    optimizer::MOI.ModelLike,
    jp::JuniperProblem;
    fix_start = false,
)
    ps = jp.options.log_levels

    jp.model = MOI.instantiate(jp.nl_solver, with_bridge_type = Float64)
    index_map = MOI.copy_to(jp.model, IntegerRelaxation(optimizer))
    # all continuous we solve relaxation first
    return jp.x = [
        index_map[vi] for vi in MOI.get(optimizer, MOI.ListOfVariableIndices())
    ]
end

function fix_primal_start!(jp::JuniperProblem)
    lb = jp.l_var
    ub = jp.u_var
    x = jp.x
    for i in 1:jp.num_var
        if jp.var_type[i] != :Cont && lb[i] <= jp.primal_start[i] <= ub[i]
            set_bounds(
                jp.model,
                jp.x[i],
                jp.primal_start[i],
                jp.primal_start[i],
            )
        else
            MOI.set(
                jp.model,
                MOI.VariablePrimalStart(),
                jp.x[i],
                jp.primal_start[i],
            )
        end
    end
end

function unfix_primal_start!(jp::JuniperProblem)
    lb = jp.l_var
    ub = jp.u_var
    x = jp.x
    for i in 1:jp.num_var
        set_bounds(jp.model, jp.x[i], lb[i], ub[i])
        MOI.set(
            jp.model,
            MOI.VariablePrimalStart(),
            jp.x[i],
            jp.primal_start[i],
        )
    end
end

function solve_root_incumbent_model(jp::JuniperProblem)
    MOI.optimize!(jp.model)
    status = MOI.get(jp.model, MOI.TerminationStatus())
    incumbent = nothing
    ps = jp.options.log_levels
    if state_is_optimal(status; allow_almost = jp.options.allow_almost_solved)
        # set incumbent
        objval = MOI.get(jp.model, MOI.ObjectiveValue())
        solution = MOI.get(jp.model, MOI.VariablePrimal(), jp.x)
        if are_type_correct(
            solution,
            jp.var_type,
            jp.disc2var_idx,
            jp.options.atol,
        )
            if only_almost_solved(status) &&
               jp.options.allow_almost_solved_integral
                @warn "Start value incumbent only almost locally solved. Disable with `allow_almost_solved_integral=false`"
            end

            incumbent = Incumbent(objval, solution, only_almost_solved(status))
            (:All in ps || :Info in ps) &&
                println("Incumbent using start values: $objval")
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
    while !state_is_optimal(
                  jp.relaxation_status;
                  allow_almost = jp.options.allow_almost_solved,
              ) &&
              restarts < max_restarts &&
              time() - jp.start_time < jp.options.time_limit

        # TODO freemodel for Knitro
        restart_values = generate_random_restart(jp)
        jp.options.debug && debug_restart_values(jp.debugDict, restart_values)
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

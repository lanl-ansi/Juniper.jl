function add_constrained_var(model, set::MOI.AbstractScalarSet)
    vi, con_idx = MOI.add_constrained_variable(model, set)
    return vi
end

function add_constrained_var(model, set::MOI.AbstractVectorSet)
    vis, con_idx = MOI.add_constrained_variables(model, set)
    return vis
end

"""
    generate_mip(optimizer, m, nlp_sol, tabu_list, start_fpump)

Generate a mip using the linear constraints of the original model
TODO: This can include quadratic constraints when the mip_solver supports them
Minimize the distance to nlp_sol and avoid using solutions inside the tabu list
"""
function generate_mip(optimizer, m, nlp_sol, tabu_list, start_fpump)
    mip_optimizer = MOI.instantiate(m.mip_solver, with_bridge_type=Float64)
    index_map = MOI.copy_to(mip_optimizer, NoObjectiveFilter(LinearFilter(optimizer)))
    x = [index_map[vi] for vi in MOI.get(optimizer, MOI.ListOfVariableIndices())]
    for i in 1:m.num_var
        if m.var_type[i] == :Bin
            if m.l_var[i] > 0
                set_bounds(mip_optimizer, x[i], 1.0, 1.0)
            elseif m.u_var[i] < 1
                set_bounds(mip_optimizer, x[i], 0.0, 0.0)
            end
        end
    end

    mabsx = add_constrained_var(mip_optimizer, MOI.Nonnegatives(m.num_disc_var))
    for (mabsxi, vi) in zip(mabsx, m.disc2var_idx)
        MOI.add_constraint(
            mip_optimizer,
            MOIU.operate(vcat, Float64, mabsxi, x[vi] - nlp_sol[vi]),
            MOI.NormOneCone(2)
        )
    end

    # How long is the tabu list
    num_sols = 0
    for i=1:tabu_list.length
        if !isnan(tabu_list.sols[i][1])
            num_sols += 1
        else
            break
        end
    end

    # If there solutions in the tabu list => avoid them
    if num_sols > 0
        z1 = [add_constrained_var(mip_optimizer, MOI.ZeroOne()) for _=1:m.num_disc_var, _=1:num_sols]
        z2 = [add_constrained_var(mip_optimizer, MOI.ZeroOne()) for _=1:m.num_disc_var, _=1:num_sols]
        v = tabu_list.sols
        for k=1:num_sols, j=1:m.num_disc_var
            i = m.disc2var_idx[j]
            lbi = m.l_var[i] > typemin(Int64) ? m.l_var[i] : typemin(Int64)
            ubi = m.u_var[i] < typemax(Int64) ? m.u_var[i] : typemax(Int64)
            MOI.add_constraint(mip_optimizer, 1.0z1[j,k] + 1.0z2[j,k], MOI.LessThan(1.0))
            MOI.add_constraint(mip_optimizer, MA.@rewrite((lbi - v[k][j]) * z1[j,k] + z2[j,k] - x[i]), MOI.LessThan(-v[k][j]))
            MOI.add_constraint(mip_optimizer, MA.@rewrite((ubi - v[k][j]) * z2[j,k] - z1[j,k] - x[i]), MOI.GreaterThan(-v[k][j]))
        end
        for k=1:num_sols
            MOI.add_constraint(mip_optimizer, sum(1.0z1) + sum(1.0z2), MOI.GreaterThan(1.0))
        end
    end

    MOI.set(mip_optimizer, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    obj = sum(1.0mabsx)
    MOI.set(mip_optimizer, MOI.ObjectiveFunction{typeof(obj)}(), obj)

    # Break the mip solver if it takes too long or throw a warning when this option isn't available
    current_time = time()-start_fpump
    time_left = m.options.feasibility_pump_time_limit-current_time
    time_left < 0 && (time_left = 1.0)

    # set time limit if supported
    old_time_limit = set_time_limit!(mip_optimizer, time_left)

    MOI.optimize!(mip_optimizer)
    status = MOI.get(mip_optimizer, MOI.TerminationStatus())

    # reset time limit
    set_time_limit!(mip_optimizer, old_time_limit)

    obj_val = NaN
    values = fill(NaN,m.num_var)
    if state_is_optimal(status; allow_almost=m.options.allow_almost_solved)
        obj_val = MOI.get(mip_optimizer, MOI.ObjectiveValue())

        # round mip values
        values = MOI.get(mip_optimizer, MOI.VariablePrimal(), x)
        for i=1:m.num_disc_var
            vi = m.disc2var_idx[i]
            values[vi] = round(values[vi])
        end
    end


    return status, values, obj_val
end

"""
    generate_nlp(optimizer, m, mip_sol, start_fpump; random_start=false)

Generates the original nlp but changes the objective to minimize the distance to the mip solution
"""
function generate_nlp(optimizer, m, mip_sol, start_fpump; random_start=false)
    nlp_optimizer = MOI.instantiate(m.nl_solver, with_bridge_type=Float64)
    index_map = MOI.copy_to(nlp_optimizer, NoObjectiveFilter(IntegerRelaxation(optimizer)))
    x = [index_map[vi] for vi in MOI.get(optimizer, MOI.ListOfVariableIndices())]
    if random_start
        restart_values = generate_random_restart(m)
    else
        restart_values = mip_sol
    end
    MOI.set(nlp_optimizer, MOI.VariablePrimalStart(), x, restart_values)

    MOI.set(nlp_optimizer, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    obj = sum((x[m.disc2var_idx[i]] - mip_sol[m.disc2var_idx[i]])^2 for i=1:m.num_disc_var)
    MOI.set(nlp_optimizer, MOI.ObjectiveFunction{typeof(obj)}(), obj)

    current_time = time()-start_fpump
    time_left = m.options.feasibility_pump_time_limit-current_time
    time_left < 0 && (time_left = 1.0)

    # set time limit if supported
    old_time_limit = set_time_limit!(nlp_optimizer, time_left)

    MOI.optimize!(nlp_optimizer)
    status = MOI.get(nlp_optimizer, MOI.TerminationStatus())

    set_time_limit!(nlp_optimizer, old_time_limit)

    nlp_obj = NaN
    nlp_sol = fill(NaN, m.num_var)
    if state_is_optimal(status; allow_almost=m.options.allow_almost_solved)
        nlp_obj = MOI.get(nlp_optimizer, MOI.ObjectiveValue())
        nlp_sol = MOI.get(nlp_optimizer, MOI.VariablePrimal(), x)
    end

    return status, nlp_sol, nlp_obj
end

"""
    generate_real_nlp(optimizer, m, sol; random_start=false)

Generate the original nlp and get the objective for that
"""
function generate_real_nlp(optimizer, m, sol; random_start=false)
    if m.num_var == m.num_disc_var
        return MOI.OPTIMAL, sol, evaluate_objective(optimizer, m, sol)
    end

    nlp_optimizer = MOI.instantiate(m.nl_solver, with_bridge_type=Float64)
    vis = collect(MOI.get(optimizer, MOI.ListOfVariableIndices()))
    disc_vals = Dict(vis[i] => sol[i] for i in m.disc2var_idx)
    index_map = MOI.copy_to(nlp_optimizer, FixVariables(IntegerRelaxation(optimizer), disc_vals))
    x = [index_map[vi] for vi in MOI.get(optimizer, MOI.ListOfVariableIndices())]

    if random_start
        restart_values = generate_random_restart(m)
        for i=1:m.num_var
            if m.var_type[i] == :Cont
                MOI.set(nlp_optimizer, MOI.VariablePrimalStart(), x[i], restart_values[i])
            else
                # discrete are fixed anyway
                MOI.set(nlp_optimizer, MOI.VariablePrimalStart(), x[i], nothing)
            end
        end
    end

    MOI.optimize!(nlp_optimizer)
    status = MOI.get(nlp_optimizer, MOI.TerminationStatus())

    obj_val = NaN
    real_sol = fill(NaN,m.num_var)
    if state_is_optimal(status; allow_almost=m.options.allow_almost_solved)
        obj_val = MOI.get(nlp_optimizer, MOI.ObjectiveValue())
        real_sol = MOI.get(nlp_optimizer, MOI.VariablePrimal(), x)
    end

    return status, real_sol, obj_val
end

"""
    add!(t::TabuList, m, sol)

Add a solution to the tabu list (includes only the discrete variables)
"""
function add!(t::TabuList, m, sol)
    t.sols[t.pointer] = [sol[m.disc2var_idx[i]] for i=1:m.num_disc_var]
    t.pointer += 1
    if t.pointer > t.length
        t.pointer = 1
    end
end

function get_fp_table(mip_obj,nlp_obj,t, fields, field_chars, catol)
    ln = ""
    i = 1
    arr = []
    for f in fields
        val = ""
        if f == :MIPobj
            if isnan(mip_obj)
                val = "-"
            else
                digits = 4
                if mip_obj < 1e-4
                    digits = convert(Int,floor(log10(1/catol)))
                end
                val = string(round(mip_obj; digits=digits))
            end
        elseif f == :NLPobj
            if isnan(nlp_obj)
                val = "-"
            else
                digits = 4
                if nlp_obj < 1e-4
                    digits = convert(Int,floor(log10(1/catol)))
                end
                val = string(round(nlp_obj; digits=digits))
            end
        elseif f == :Time
            val = string(round(t; digits=1))
        end

        if length(val) > field_chars[i]
            # too long to display shouldn't happen normally but is better than error
            # if it happens
            val = "t.l."
        end

        padding = field_chars[i]-length(val)
        ln *= repeat(" ",trunc(Int, floor(padding/2)))
        ln *= val
        ln *= repeat(" ",trunc(Int, ceil(padding/2)))
        push!(arr,val)
        i += 1
    end
    return ln, arr
end

"""
    fpump(optimizer, m)

Run the feasibility pump
"""
function fpump(optimizer, m)

    Random.seed!(JUNIPER_RNG, m.options.seed)

    if are_type_correct(m.relaxation_solution, m.var_type, m.disc2var_idx, m.options.atol)
        return Incumbent(m.relaxation_objval, m.relaxation_solution, only_almost_solved(m.status))
    end

    start_fpump = time()
    nlp_sol = m.relaxation_solution
    nlp_obj = 1 # should be not 0 for while
    c = 1
    tabu_list = TabuList()
    mip_sols = Dict{UInt64,Bool}()
    tabu_list.length = m.options.tabu_list_length
    tabu_list.pointer = 1
    tabu_list.sols = []
    for i=1:tabu_list.length
        push!(tabu_list.sols, NaN*ones(m.num_disc_var))
    end

    last_table_arr = []
    fields = []
    field_chars = []
    ps = m.options.log_levels
    # Print table init
    if check_print(ps,[:Table])
        fields, field_chars = [:MIPobj,:NLPobj,:Time], [20,20,5]
        print_table_header(fields,field_chars)
    end

    real_status = MOI.OPTIMIZE_NOT_CALLED
    fix = false
    nlp_status = :Error
    iscorrect = false
    tl = m.options.feasibility_pump_time_limit
    # the tolerance can be changed => current atol
    catol = m.options.atol
    atol_counter = 0
    while !are_type_correct(nlp_sol, m.var_type, m.disc2var_idx, catol) && time()-start_fpump < tl &&
        time()-m.start_time < m.options.time_limit

        # generate a mip or just round if no linear constraints
        if any(FS -> FS[1] != MOI.VariableIndex, MOI.get(LinearFilter(optimizer), MOI.ListOfConstraintTypesPresent()))
            mip_status, mip_sol, mip_obj = generate_mip(optimizer, m, nlp_sol, tabu_list, start_fpump)
        else
            # if no linear constraints just round the discrete variables
            mip_obj = NaN
            mip_sol = copy(nlp_sol)
            mip_status = MOI.OPTIMAL
            for vi=1:m.num_disc_var
                vidx = m.disc2var_idx[vi]
                mip_sol[vidx] = round(mip_sol[vidx])
            end
        end
        if mip_status != MOI.OPTIMAL
            @warn "MIP couldn't be solved to optimality. Terminated with status: "*string(mip_status)
            break
        end

        # If a cycle is detected which wasn't able to prevent by the tabu list (maybe too short)
        if haskey(mip_sols, hash(mip_sol))
            @warn "Cycle detected"
            break
        end
        add!(tabu_list, m, mip_sol)
        mip_sols[hash(mip_sol)] = true

        nlp_status, nlp_sol, nlp_obj = generate_nlp(optimizer, m, mip_sol, start_fpump)
        if !state_is_optimal(nlp_status; allow_almost=m.options.allow_almost_solved)
            cnlpinf = 0
            while cnlpinf < m.options.num_resolve_nlp_feasibility_pump && !state_is_optimal(nlp_status; allow_almost=m.options.allow_almost_solved) &&
                time()-start_fpump < tl && time()-m.start_time < m.options.time_limit
                nlp_status, nlp_sol, nlp_obj = generate_nlp(optimizer, m, mip_sol, start_fpump; random_start=true)
                cnlpinf += 1
            end
            if !state_is_optimal(nlp_status; allow_almost=m.options.allow_almost_solved)
                @warn "NLP couldn't be solved to optimality"
                if check_print(ps,[:Table])
                    print_fp_table(mip_obj, NaN, time()-start_fpump, fields, field_chars, catol)
                end
                break
            end
        end

        if check_print(ps,[:Table])
            print_fp_table(mip_obj, nlp_obj, time()-start_fpump, fields, field_chars, catol)
        end

        # if the current tolerance was nearly reached 5 times
        # => If reasonable should be an option
        if atol_counter >= m.options.feasibility_pump_tolerance_counter
            catol *= 10
            @warn "FPump tolerance changed to: ",catol
            atol_counter = 0
        end

        # if the difference is near 0 => try to improve the obj by using the original obj
        # set atol for type correct to a low value as it is checked with real_nlp anyway
        if are_type_correct(nlp_sol, m.var_type, m.disc2var_idx, catol*1000) || isapprox(nlp_obj, 0.0; atol=catol)
            real_status,real_sol, real_obj = generate_real_nlp(optimizer, m, mip_sol)
            cnlpinf = 0
            while cnlpinf < m.options.num_resolve_nlp_feasibility_pump && !state_is_optimal(real_status; allow_almost=m.options.allow_almost_solved) &&
                time()-start_fpump < tl && time()-m.start_time < m.options.time_limit
                real_status,real_sol, real_obj = generate_real_nlp(optimizer, m, mip_sol; random_start=true)
                cnlpinf += 1
            end
            if state_is_optimal(real_status) || (only_almost_solved(real_status) && m.options.allow_almost_solved_integral)
                if only_almost_solved(real_status)
                    @warn "Integral feasible point only almost solved. Disable with `allow_almost_solved_integral=false`"
                end
                nlp_obj = real_obj
                nlp_sol = real_sol
                iscorrect = true
                break
            elseif are_type_correct(nlp_sol, m.var_type, m.disc2var_idx, catol)
                real_status = MOI.LOCALLY_SOLVED
                nlp_obj = evaluate_objective(optimizer, m, nlp_sol)
                iscorrect = true
                @warn "Real objective wasn't solved to optimality"
                break
            end
        end
        if !isapprox(nlp_obj, 0.0; atol=catol) && isapprox(nlp_obj, 0.0; atol=10*catol)
            atol_counter += 1
        else
            atol_counter = 0
        end
        c += 1
    end

    if check_print(ps,[:Table])
        println()
    end

    if check_print(ps,[:Info])
        println("FP: ", time()-start_fpump, " s")
        println("FP: ", c == 1 ? "$c round" : "$c rounds")
    end
    m.fpump_info = Dict{Symbol,Float64}()
    m.fpump_info[:time] = time()-start_fpump
    m.fpump_info[:rounds] = c

    if iscorrect
        check_print(ps,[:Info]) && println("FP: Obj: ", nlp_obj)
        m.fpump_info[:obj] = nlp_obj
        m.fpump_info[:gap] = abs(m.relaxation_objval-nlp_obj)/abs(nlp_obj)
        return Incumbent(nlp_obj, nlp_sol, only_almost_solved(real_status))
    end

    m.fpump_info[:obj] = NaN
    m.fpump_info[:gap] = NaN
    check_print(ps,[:Info]) && println("FP: No integral solution found")

    return nothing
end

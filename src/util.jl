#=
    Used from https://github.com/lanl-ansi/Alpine.jl
=# 
function expr_dereferencing!(expr, m)
    for i in 2:length(expr.args)
        if isa(expr.args[i], Union{Float64,Int64})
            k = 0
        elseif expr.args[i].head == :ref
            @assert isa(expr.args[i].args[2], MOI.VariableIndex)
            expr.args[i] = JuMP.VariableRef(m, expr.args[i].args[2])
        elseif expr.args[i].head == :call
            expr_dereferencing!(expr.args[i], m)
        else
            error("expr_dereferencing :: Unexpected term in expression tree.")
        end
    end
end

function generate_random_restart(m; cont=true)
    values = []
    for i=1:m.num_var
        lbi_def = true
        ubi_def = true
        if m.l_var[i] > typemin(Int64) 
            lbi = m.l_var[i]
        else
            lbi = typemin(Int64)
            lbi_def = false
        end

        if m.u_var[i] < typemax(Int64) 
            ubi = m.u_var[i]
        else
            ubi = typemin(Int64)
            ubi_def = false
        end

        if !ubi_def && !lbi_def 
            ubi = 10
            lbi = -10
        elseif !ubi_def
            ubi = lbi+20
        elseif !lbi_def
            lbi = ubi-20
        end             

        if m.var_type[i] == :Cont || cont
            push!(values,(ubi-lbi)*rand()+lbi)
        else
            lbi = Int(round(lbi))
            ubi = Int(round(ubi))
            push!(values, rand(lbi:ubi))
        end
    end
    return values
end

"""
    get_reasonable_disc_vars(node, var_type, disc_vars, disc2var_idx, atol)

Get all discrete variables which aren't close to discrete yet based on atol 
"""
function get_reasonable_disc_vars(node, var_type, disc_vars, disc2var_idx, atol)
    reasonable_disc_vars = zeros(Int64,0)
    for i=1:disc_vars
        idx = disc2var_idx[i]
        u_b = node.u_var[idx]
        l_b = node.l_var[idx]
        if isapprox(u_b,l_b; atol=atol) || is_type_correct(node.solution[idx],var_type[idx],atol)
            continue
        end
        push!(reasonable_disc_vars,i)
    end
    return reasonable_disc_vars
end

function get_type_dict(obj)
    T = typeof(obj)
    type_dict = Dict{Symbol,Type}()
    for (name, typ) in zip(fieldnames(T), T.types)
        type_dict[name] = typ
    end
    return type_dict
end

"""
    is_global_status(state::MOI.TerminationStatusCode)

Returns true if either ALMOST_OPTIMAL, OPTIMAL or INFEASIBLE and false otherwise
"""
function is_global_status(state::MOI.TerminationStatusCode)
    return state == MOI.ALMOST_OPTIMAL || state == MOI.OPTIMAL || state == MOI.INFEASIBLE
end

"""
    only_almost_solved(state::MOI.TerminationStatusCode)

Returns true if either ALMOST_OPTIMAL or ALMOST_LOCALLY_SOLVED
"""
function only_almost_solved(state::MOI.TerminationStatusCode)
    return state == MOI.ALMOST_OPTIMAL || state == MOI.ALMOST_LOCALLY_SOLVED
end

"""
    state_is_optimal(state::MOI.TerminationStatusCode; allow_almost=false)

Returns true if either optimal or locally solved. If allow_almost then check for `ALMOST_LOCALLY_SOLVED` and `ALMOST_OPTIMAL`
"""
function state_is_optimal(state::MOI.TerminationStatusCode; allow_almost=false)
    return state == MOI.OPTIMAL || state == MOI.LOCALLY_SOLVED || 
            (allow_almost && state == MOI.ALMOST_LOCALLY_SOLVED) || (allow_almost && state == MOI.ALMOST_OPTIMAL)
end

"""
    state_is_infeasible(state::MOI.TerminationStatusCode)

Returns true if either infeasible or locally infeasible
"""
function state_is_infeasible(state::MOI.TerminationStatusCode)
    return state == MOI.INFEASIBLE || state == MOI.LOCALLY_INFEASIBLE
end

"""
    add_obj_constraint(jp::JuniperProblem, rhs::Float64)

Add a constraint for the objective based on whether the objective is linear/quadratic or non linear.
If the objective sense is :MIN than add objective <= rhs else objective >= rhs
"""
function add_obj_constraint(jp::JuniperProblem, rhs::Float64)
    if jp.has_nl_objective
        obj_expr = MOI.objective_expr(jp.nlp_evaluator)
        if jp.obj_sense == :Min
            obj_constr = Expr(:call, :<=, obj_expr, rhs)
        else
            obj_constr = Expr(:call, :>=, obj_expr, rhs)
        end
        Juniper.expr_dereferencing!(obj_constr, jp.model)
        JuMP.add_NL_constraint(jp.model, obj_constr)
    else # linear or quadratic
        backend = JuMP.backend(jp.model);
        if isa(jp.objective, MOI.SingleVariable)
            if jp.obj_sense == :Min
                JuMP.set_upper_bound(jp.x[jp.objective.variable.value], rhs)
            else
                JuMP.set_lower_bound(jp.x[jp.objective.variable.value], rhs)
            end
        else
            if jp.obj_sense == :Min
                MOI.add_constraint(backend, jp.objective, MOI.LessThan(rhs))
            else
                MOI.add_constraint(backend, jp.objective, MOI.GreaterThan(rhs))
            end
        end
    end
end


"""
    evaluate_objective(optimizer::MOI.AbstractOptimizer, jp::JuniperProblem, xs::Vector{Float64})

If no objective exists => return 0
Evaluate the objective whether it is non linear, linear or quadratic
"""
function evaluate_objective(optimizer::MOI.AbstractOptimizer, jp::JuniperProblem, xs::Vector{Float64})
    if optimizer.nlp_data.has_objective
        return MOI.eval_objective(optimizer.nlp_data.evaluator, xs)
    elseif optimizer.objective !== nothing
        return MOIU.eval_variables(vi -> xs[vi.value], optimizer.objective)
    else 
        return 0.0
    end
end

"""
    set_subsolver_option!(jp::JuniperProblem, model::JuMP.model, type_of_subsolver::String,
                          subsolver_name::String, param::Symbol, change::Pair)

Set the optimizer of the model if the subsolver_name is part of the name of the optimizer.
Change the option `param` of the sub_solver. i.e `change=0.1 => 1e-5` means the default 
normally is 0.1 and will be changed to 1e-5.
Normally reset_subsolver_option! should be run to return to reset this change after `optimize!`
Return the previous value of the param option or if not previously set return change.first   
"""
function set_subsolver_option!(jp::JuniperProblem, model::JuMP.Model, type_of_subsolver::String,
                              subsolver_name::String, param::Symbol, change::Pair)

    old_value = change.first
    if type_of_subsolver == "nl"
        sub_solver = getfield(jp, :nl_solver)
        sub_solver_options = getfield(jp, :nl_solver_options)
    elseif type_of_subsolver == "mip"
        sub_solver = getfield(jp, :mip_solver)
        sub_solver_options = getfield(jp, :mip_solver_options)
    end

    if occursin(subsolver_name, string(sub_solver))
        overwritten = false
        for i=1:length(sub_solver_options)
            if sub_solver_options[i][1] == param
                old_value = sub_solver_options[i][2]
                sub_solver_options[i] = (param, change.second)
                overwritten = true
                break
            end
        end
        if !overwritten
            push!(sub_solver_options, (param, change.second))
        end
        sub_solver = with_optimizer(sub_solver.constructor; sub_solver_options...)
    end

    if type_of_subsolver == "nl"
        setfield!(jp, :nl_solver, sub_solver)
        setfield!(jp, :nl_solver_options, sub_solver_options)
        JuMP.set_optimizer(model, jp.nl_solver)  
    elseif type_of_subsolver == "mip"
        setfield!(jp, :mip_solver, sub_solver)
        setfield!(jp, :mip_solver_options, sub_solver_options)
        JuMP.set_optimizer(model, jp.mip_solver) 
    end

    
    return old_value
end

"""
    reset_subsolver_option!(jp::JuniperProblem, type_of_subsolver::String,
                            subsolver_name::String, param::Symbol, value)

Resets the subsolver option `param` to `value` if the subsolver for the type i.e "nl" matches
`subsolver_name`. `value` is normally get by calling `set_subsolver_option!`
"""
function reset_subsolver_option!(jp::JuniperProblem, type_of_subsolver::String,
                                subsolver_name::String, param::Symbol, value)
    if type_of_subsolver == "nl"
        sub_solver = getfield(jp, :nl_solver)
        sub_solver_options = getfield(jp, :nl_solver_options)
    elseif type_of_subsolver == "mip"
        sub_solver = getfield(jp, :mip_solver)
        sub_solver_options = getfield(jp, :mip_solver_options)
    end

    if occursin(subsolver_name, string(sub_solver))
        for i=1:length(sub_solver_options)
            if sub_solver_options[i][1] == param
                sub_solver_options[i] = (param, value)
                break
            end
        end
        sub_solver = with_optimizer(sub_solver.constructor; sub_solver_options...)
    end

    if type_of_subsolver == "nl"
        setfield!(jp, :nl_solver, sub_solver)
        setfield!(jp, :nl_solver_options, sub_solver_options)
    elseif type_of_subsolver == "mip"
        setfield!(jp, :mip_solver, sub_solver)
        setfield!(jp, :mip_solver_options, sub_solver_options)
    end
end

"""
    optimize_get_status_backend(model::JuMP.Model; solver::Union{Nothing,JuMP.OptimizerFactory}=nothing) 

Run optimize! and get the status and the backend
"""
function optimize_get_status_backend(model::JuMP.Model; solver::Union{Nothing,JuMP.OptimizerFactory}=nothing) 
    if solver == nothing
        JuMP.optimize!(model)
    else
        JuMP.optimize!(model, solver)
    end
    backend = JuMP.backend(model)
    status = MOI.get(backend, MOI.TerminationStatus()) 
    return status, backend
end
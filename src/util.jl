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
            push!(values,(ubi-lbi)*rand(JUNIPER_RNG;)+lbi)
        else
            lbi = Int(round(lbi))
            ubi = Int(round(ubi))
            push!(values, rand(JUNIPER_RNG,lbi:ubi))
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
        if isa(jp.objective, MOI.SingleVariable)
            i = findfirst(isequal(jp.objective.variable), jp.x)
            set = MOI.get(jp.model, MOI.ConstraintSet(), jp.cx[i])
            if jp.obj_sense == :Min
                set = MOI.Interval(set.lower, rhs)
            else
                set = MOI.Interval(rhs, set.upper)
            end
            set = MOI.get(jp.model, MOI.ConstraintSet(), jp.cx[i], set)
        else
            if jp.obj_sense == :Min
                MOI.add_constraint(jp.model, jp.objective, MOI.LessThan(rhs))
            else
                MOI.add_constraint(jp.model, jp.objective, MOI.GreaterThan(rhs))
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
                          subsolver_name::String, param::String, change::Pair)

Set the optimizer of the model if the subsolver_name is part of the name of the optimizer.
Change the option `param` of the sub_solver. i.e `change=0.1 => 1e-5` means the default
normally is 0.1 and will be changed to 1e-5.
Normally reset_subsolver_option! should be run to return to reset this change after `optimize!`
Return the previous value of the param option or if not previously set return change.first
"""
function set_subsolver_option!(model::MOI.ModelLike, attr::MOI.AbstractOptimizerAttribute, change::Pair)
    prev = try
        MOI.get(model, attr)
    catch
        change.first
    end
    MOI.set(model, attr, change.second)
    return prev
end
function set_subsolver_option!(model::MOI.ModelLike, attr::MOI.AbstractOptimizerAttribute, change)
    MOI.set(model, attr, change)
    return
end
function set_subsolver_option!(model::MOI.ModelLike, subsolver_name::String, attr::MOI.AbstractOptimizerAttribute, change)
    if occursin(subsolver_name, MOI.get(model, MOI.SolverName()))
        return set_subsolver_option!(model, attr, change)
    end
    return value
end
function set_subsolver_option!(model::MOI.ModelLike, subsolver_name::String, param::String, change)
    set_subsolver_option!(model, subsolver_name, MOI.RawParameter(param), change)
end

function set_time_limit!(optimizer, time_limit::Union{Nothing,Float64})
    old_time_limit = Inf
    if MOI.supports(optimizer, MOI.TimeLimitSec())
        old_time_limit = MOI.get(optimizer, MOI.TimeLimitSec())
        MOI.set(optimizer, MOI.TimeLimitSec(), time_limit)
    end
    return old_time_limit
end

# TODO maybe move to MOI.Utilities alongside `MOIU.get_bounds`
function set_bounds(model::MOI.ModelLike, vi::MOI.VariableIndex, lower::T, upper::T) where T
    xval = vi.value
    c_lt = MOI.ConstraintIndex{MOI.SingleVariable,MOI.LessThan{T}}(xval)
    c_gt = MOI.ConstraintIndex{MOI.SingleVariable,MOI.GreaterThan{T}}(xval)
    c_int = MOI.ConstraintIndex{MOI.SingleVariable,MOI.Interval{T}}(xval)
    c_eq = MOI.ConstraintIndex{MOI.SingleVariable,MOI.EqualTo{T}}(xval)
    if MOI.is_valid(model, c_int)
        MOI.set(model, MOI.ConstraintSet(), c_int, MOI.Interval(lower, upper))
        # It is assumed that none of the other ConstraintIndexs are valid
        return
    end
    if MOI.is_valid(model, c_eq)
        # Maybe delete it and add an interval constraint if this case is useful
        lower == upper || error("Cannot set different bounds $lower and $upper to a variables with `EqualTo` constraint.")
        MOI.set(model, MOI.ConstraintSet(), c_eq, MOI.EqualTo(lower))
        # It is assumed that none of the other `ConstraintIndex`s are valid
        return
    end
    lt = MOI.LessThan(upper)
    if MOI.is_valid(model, c_lt)
        if upper == typemin(upper)
            MOI.delete(model, c_lt)
        else
            MOI.set(model, MOI.ConstraintSet(), c_lt, lt)
        end
    elseif upper != typemax(upper)
        MOI.add_constraint(model, MOI.SingleVariable(vi), lt)
    end
    gt = MOI.GreaterThan(lower)
    if MOI.is_valid(model, c_gt)
        if lower == typemin(lower)
            MOI.delete(model, c_gt)
        else
            MOI.set(model, MOI.ConstraintSet(), c_gt, gt)
        end
    elseif lower != typemin(lower)
        MOI.add_constraint(model, MOI.SingleVariable(vi), gt)
    end
    return
end
function set_lower_bound(model::MOI.ModelLike, vi::MOI.VariableIndex, lower::T) where T
    xval = vi.value
    c_gt = MOI.ConstraintIndex{MOI.SingleVariable,MOI.GreaterThan{T}}(xval)
    c_int = MOI.ConstraintIndex{MOI.SingleVariable,MOI.Interval{T}}(xval)
    c_eq = MOI.ConstraintIndex{MOI.SingleVariable,MOI.EqualTo{T}}(xval)
    if MOI.is_valid(model, c_int)
        set = MOI.get(model, MOI.ConstraintSet(), c_int)
        MOI.set(model, MOI.ConstraintSet(), c_int, MOI.Interval(lower, set.upper))
        # It is assumed that none of the other ConstraintIndexs are valid
        return
    end
    if MOI.is_valid(model, c_eq)
        error("Cannot set lower bound to fixed variable.")
    end
    gt = MOI.GreaterThan(lower)
    if MOI.is_valid(model, c_gt)
        MOI.set(model, MOI.ConstraintSet(), c_gt, gt)
    else
        MOI.add_constraint(model, MOI.SingleVariable(vi), gt)
    end
    return
end
function set_upper_bound(model::MOI.ModelLike, vi::MOI.VariableIndex, upper::T) where T
    xval = vi.value
    c_lt = MOI.ConstraintIndex{MOI.SingleVariable,MOI.LessThan{T}}(xval)
    c_int = MOI.ConstraintIndex{MOI.SingleVariable,MOI.Interval{T}}(xval)
    c_eq = MOI.ConstraintIndex{MOI.SingleVariable,MOI.EqualTo{T}}(xval)
    if MOI.is_valid(model, c_int)
        set = MOI.get(model, MOI.ConstraintSet(), c_int)
        MOI.set(model, MOI.ConstraintSet(), c_int, MOI.Interval(set.lower, upper))
        # It is assumed that none of the other ConstraintIndexs are valid
        return
    end
    if MOI.is_valid(model, c_eq)
        error("Cannot set lower bound to fixed variable.")
    end
    lt = MOI.LessThan(upper)
    if MOI.is_valid(model, c_lt)
        MOI.set(model, MOI.ConstraintSet(), c_lt, lt)
    else
        MOI.add_constraint(model, MOI.SingleVariable(vi), lt)
    end
    return
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

struct ObjectiveConstraint{E<:MOI.AbstractNLPEvaluator} <: MOI.AbstractNLPEvaluator
    evaluator::E
    num_variables::Int
end
MOI.features_available(oc::ObjectiveConstraint) = MOI.features_available(oc.evaluator)
MOI.initialize(oc::ObjectiveConstraint, features) = MOI.initialize(oc.evaluator, features)
MOI.eval_objective(oc::ObjectiveConstraint, x) = MOI.eval_objective(oc.evaluator, x)
MOI.eval_objective_gradient(oc::ObjectiveConstraint, g, x) = MOI.eval_objective_gradient(oc.evaluator, g, c)
MOI.hessian_lagrangian_structure(oc::ObjectiveConstraint) = MOI.hessian_lagrangian_structure(oc.evaluator)
MOI.eval_hessian_lagrangian(oc::ObjectiveConstraint, H, x, σ, μ) = MOI.eval_hessian_lagrangian(oc.evaluator, H, x, σ, μ)
function MOI.eval_constraint(oc::ObjectiveConstraint, g, x)
    g[1] = MOI.eval_objective(oc.evaluator, x)
    MOI.eval_constraint(oc.evaluator, view(g, 2:length(g)), x)
    return
end
function MOI.jacobian_structure(oc::ObjectiveConstraint)
    sparsity = [(1, i) for i in 1:oc.num_variables]
    append!(sparsity, MOI.jacobian_structure(oc.evaluator))
    return sparsity
end
function MOI.eval_constraint_jacobian(oc::ObjectiveConstraint, J, x)
    MOI.eval_objective_gradient(oc.evaluator, view(J, 1:oc.num_variables), x)
    MOI.eval_constraint_jacobian(oc.evaluator, view(J, (oc.num_variables + 1):length(J)), x)
    return
end

"""
    add_obj_constraint(jp::JuniperProblem, rhs::Float64)

Add a constraint for the objective based on whether the objective is linear/quadratic or non linear.
If the objective sense is :MIN than add objective <= rhs else objective >= rhs
"""
function add_obj_constraint(jp::JuniperProblem, rhs::Float64)
    if jp.nlp_data.has_objective
        if jp.obj_sense == :Min
            lb, ub = -Inf, rhs
        else
            lb, ub = rhs, Inf
        end
        MOI.set(jp.model, MOI.NLPBlock(), MOI.NLPBlockData(
            [jp.nlp_data.constraint_bounds; MOI.NLPBoundsPair(lb, ub)],
            ObjectiveConstraint(jp.nlp_data.evaluator, MOI.get(jp.model, MOI.NumberOfVariables())),
            jp.nlp_data.has_objective,
        ))
    else # linear or quadratic
        if isa(jp.objective, MOI.SingleVariable)
            i = findfirst(isequal(jp.objective.variable), jp.x)
            vi = jp.objective.variable
            if jp.obj_sense == :Min
                set_upper_bound(jp.model, vi, rhs)
            else
                set_lower_bound(jp.model, vi, rhs)
            end
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
    # NLP objective has priority if both are set
    if jp.nlp_data !== nothing && jp.nlp_data.has_objective
        return MOI.eval_objective(jp.nlp_data.evaluator, xs)
    elseif jp.objective !== nothing
        return MOIU.eval_variables(vi -> xs[vi.value], jp.objective)
    else
        return 0.0
    end
end

"""
    set_subsolver_option!(model::MOI.ModelLike, attr::MOI.AbstractOptimizerAttribute, change::Pair)

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
end
function set_subsolver_option!(model::MOI.ModelLike, subsolver_name::String, param::String, change)
    set_subsolver_option!(model, subsolver_name, MOI.RawOptimizerAttribute(param), change)
end

function set_time_limit!(optimizer, time_limit::Union{Nothing,Float64})
    old_time_limit = Inf
    if MOI.supports(optimizer, MOI.TimeLimitSec())
        old_time_limit = MOI.get(optimizer, MOI.TimeLimitSec())
        MOI.set(optimizer, MOI.TimeLimitSec(), time_limit)
    end
    return old_time_limit
end

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
    get_reasonable_int_vars(node, var_type, int_vars, disc2var_idx, atol)

Get all discrete variables which aren't close to discrete yet based on atol 
"""
function get_reasonable_int_vars(node, var_type, int_vars, disc2var_idx, atol)
    reasonable_int_vars = zeros(Int64,0)
    for i=1:int_vars
        idx = disc2var_idx[i]
        u_b = node.u_var[idx]
        l_b = node.l_var[idx]
        if isapprox(u_b,l_b; atol=atol) || is_type_correct(node.solution[idx],var_type[idx],atol)
            continue
        end
        push!(reasonable_int_vars,i)
    end
    return reasonable_int_vars
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
    state_is_optimal(state::MOI.TerminationStatusCode)

Returns true if either optimal or locally solved
"""
function state_is_optimal(state::MOI.TerminationStatusCode)
    return state == MOI.OPTIMAL || state == MOI.LOCALLY_SOLVED
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

Evaluate the objective whether it is non linear, linear or quadratic
"""
function evaluate_objective(optimizer::MOI.AbstractOptimizer, jp::JuniperProblem, xs::Vector{Float64})
    if optimizer.nlp_data.has_objective
        return MOI.eval_objective(optimizer.nlp_data.evaluator, xs)
    else
        return MOIU.evalvariables(vi -> xs[vi.value], optimizer.objective)
    end
end
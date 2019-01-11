# TODO(odow):
# - remove the `try-catch`es once the following PR is merged
#   https://github.com/JuliaOpt/MathOptInterface.jl/pull/623

using MathOptInterface
const MOI = MathOptInterface
const MOIU = MOI.Utilities

import MathProgBase
const MPB = MathProgBase.SolverInterface

MOIU.@model(InnerModel,
    (MOI.ZeroOne, MOI.Integer),
    (MOI.EqualTo, MOI.GreaterThan, MOI.LessThan, MOI.Interval),
    (),
    (),
    (MOI.SingleVariable,),
    (MOI.ScalarAffineFunction, MOI.ScalarQuadraticFunction),
    (),
    ()
)

const MOI_SCALAR_SETS = (
    MOI.EqualTo{Float64}, MOI.GreaterThan{Float64}, MOI.LessThan{Float64},
    MOI.Interval{Float64}
)

const MOI_SCALAR_FUNCTIONS = (
    MOI.ScalarAffineFunction{Float64}, MOI.ScalarQuadraticFunction{Float64}
)

# Make `Model` a constant for use in the rest of the file.
const Model = MOIU.UniversalFallback{InnerModel{Float64}}

# Only support the constraint types defined by `InnerModel`.
function MOI.supports_constraint(model::Model, F::Type{<:MOI.AbstractFunction},
        S::Type{<:MOI.AbstractSet})
    return MOI.supports_constraint(model.model, F, S)
end

"Attribute for the MathProgBase solver."
struct MPBSolver <: MOI.AbstractOptimizerAttribute end

"Attribute for the MathProgBase status."
struct MPBSolutionAttribute <: MOI.AbstractModelAttribute end

"Struct to contain the MPB solution."
struct MPBSolution
    status::Symbol
    objective_value::Float64
    primal_solution::Dict{MOI.VariableIndex, Float64}
end

"""
    Optimizer()
"""
function Optimizer()
    model = MOIU.UniversalFallback(InnerModel{Float64}())
    MOI.set(model, MPBSolver(),
        # AmplNLSolver(solver_command, options, filename = filename)
    )
    return model
end

Base.show(io::IO, ::Model) = println(io, "A Juniper model")

MOI.get(::Model, ::MOI.SolverName) = "Juniper"

# We re-define is_empty and empty! to prevent the universal fallback from
# deleting the solver that we are caching in it.
function MOI.is_empty(model::Model)
    return MOI.is_empty(model.model) &&
        isempty(model.constraints) &&
        isempty(model.modattr) &&
        isempty(model.varattr) &&
        isempty(model.conattr) &&
        length(model.optattr) == 1 &&
        haskey(model.optattr, MPBSolver())
end

function MOI.empty!(model::Model)
    MOI.empty!(model.model)
    empty!(model.constraints)
    model.nextconstraintid = 0
    empty!(model.con_to_name)
    model.name_to_con = nothing
    empty!(model.modattr)
    empty!(model.varattr)
    empty!(model.conattr)
    mpb_solver = model.optattr[MPBSolver()]
    empty!(model.optattr)
    model.optattr[MPBSolver()] = mpb_solver
    return
end

set_to_bounds(set::MOI.LessThan) = (-Inf, set.upper)
set_to_bounds(set::MOI.GreaterThan) = (set.lower, Inf)
set_to_bounds(set::MOI.EqualTo) = (set.value, set.value)
set_to_bounds(set::MOI.Interval) = (set.lower, set.upper)
set_to_cat(set::MOI.ZeroOne) = :Bin
set_to_cat(set::MOI.Integer) = :Int

struct NLPEvaluator{T} <: MPB.AbstractNLPEvaluator
    inner::T
    variable_map::Dict{MOI.VariableIndex, Int}
    num_inner_con::Int
    objective_expr::Union{Nothing, Expr}
    scalar_constraint_expr::Vector{Expr}
end

"""
MathProgBase expects expressions with variables denoted by `x[i]` for contiguous
`i`. However, JuMP 0.19 creates expressions with `x[MOI.VariableIndex(i)]`. So
we have to recursively walk the expression replacing instances of
MOI.VariableIndex by a corresponding integer.
"""
function replace_variableindex_by_int(variable_map, expr::Expr)
    for (i, arg) in enumerate(expr.args)
        expr.args[i] = replace_variableindex_by_int(variable_map, arg)
    end
    return expr
end
function replace_variableindex_by_int(variable_map, expr::MOI.VariableIndex)
    return variable_map[expr]
end
replace_variableindex_by_int(variable_map, expr) = expr

# ==============================================================================
# In the next section we match up the MPB nonlinear functions to the MOI API.
# This is basically just a rename.
function MPB.initialize(d::NLPEvaluator, requested_features::Vector{Symbol})
    if d.inner !== nothing
        MOI.initialize(d.inner, requested_features)
    end
    return
end

function MPB.features_available(d::NLPEvaluator)
    if d.inner !== nothing
        return MOI.features_available(d.inner)
    else
        return [:ExprGraph]
    end
end

function MPB.obj_expr(d::NLPEvaluator)
    # d.objective_expr is a SingleVariable, ScalarAffineFunction, or a
    # ScalarQuadraticFunction from MOI. If it is unset, it will be `nothing` (we
    # enforce this when creating the NLPEvaluator in `optimize!`).
    if d.objective_expr !== nothing
        return d.objective_expr
    elseif d.inner !== nothing
        # If d.objective_expr === nothing, then the objective must be nonlinear.
        # Query it from the inner NLP evaluator.
        expr = MOI.objective_expr(d.inner)
        return replace_variableindex_by_int(d.variable_map, expr)
    else
        return :(0.0)
    end
end

function MPB.constr_expr(d::NLPEvaluator, i)
    # There are two types of constraints in the model:
    # i = 1, 2, ..., d.num_inner_con are the nonlinear constraints. We access
    # these from the inner NLP evaluator.
    # i = d.num_inner_con + 1, d.num_inner_con + 2, ..., N are the linear or
    # quadratic constraints added by MOI.
    if i <= d.num_inner_con
        expr = MOI.constraint_expr(d.inner, i)
        return replace_variableindex_by_int(d.variable_map, expr)
    else
        return d.scalar_constraint_expr[i - d.num_inner_con]
    end
end

# ==============================================================================
# In the next section, we need to convert MOI functions and sets into expression
# graphs. First, we convert functions (SingleVariable, ScalarAffine, and
# ScalarQuadratic) into expression graphs.
function func_to_expr_graph(func::MOI.SingleVariable, variable_map)
    return Expr(:ref, :x, variable_map[func.variable])
end

function func_to_expr_graph(func::MOI.ScalarAffineFunction, variable_map)
    expr = Expr(:call, :+, func.constant)
    for term in func.terms
        push!(expr.args, Expr(:call, :*, term.coefficient,
            Expr(:ref, :x, variable_map[term.variable_index])
        ))
    end
    return expr
end

function func_to_expr_graph(func::MOI.ScalarQuadraticFunction, variable_map)
    expr = Expr(:call, :+, func.constant)
    for term in func.affine_terms
        push!(expr.args, Expr(:call, :*, term.coefficient,
            Expr(:ref, :x, variable_map[term.variable_index])
        ))
    end
    for term in func.quadratic_terms
        index_1 = variable_map[term.variable_index_1]
        index_2 = variable_map[term.variable_index_2]
        coef = term.coefficient
        # MOI defines quadratic as 1/2 x' Q x :(
        if index_1 == index_2
            coef *= 0.5
        end
        push!(expr.args, Expr(:call, :*, coef,
            Expr(:ref, :x, index_1),
            Expr(:ref, :x, index_2)
        ))
    end
    return expr
end

# Next, we need to combine a function `Expr` with a MOI set into a comparison.

function funcset_to_expr_graph(func::Expr, set::MOI.LessThan)
    return Expr(:call, :<=, func, set.upper)
end

function funcset_to_expr_graph(func::Expr, set::MOI.GreaterThan)
    return Expr(:call, :>=, func, set.lower)
end

function funcset_to_expr_graph(func::Expr, set::MOI.EqualTo)
    return Expr(:call, :(==), func, set.value)
end

function funcset_to_expr_graph(func::Expr, set::MOI.Interval)
    return Expr(:comparison, set.lower, :<=, func, :<=, set.upper)
end

# A helper function that converts a MOI function and MOI set into a comparison
# expression.

function moi_to_expr_graph(func, set, variable_map)
    func_expr = func_to_expr_graph(func, variable_map)
    return funcset_to_expr_graph(func_expr, set)
end

# ==============================================================================
# Now we're ready to optimize model.
function MOI.optimize!(model::Model)
    # Extract the MathProgBase solver from the model. Recall it's a MOI
    # attribute which we're careful not to delete on calls to `empty!`.
    mpb_solver = MOI.get(model, MPBSolver())

    # Get the optimzation sense.
    opt_sense = MOI.get(model, MOI.ObjectiveSense())
    sense = opt_sense == MOI.MAX_SENSE ? :Max : :Min

    # Get the NLPBlock from the model.
    nlp_block = try
        MOI.get(model, MOI.NLPBlock())
    catch
        nothing
    end

    # ==========================================================================
    # Extract the nonlinear constraint bounds. We're going to append to these
    # g_l and g_u vectors later.
    num_con = 0
    moi_nlp_evaluator = nlp_block !== nothing ? nlp_block.evaluator : nothing
    if nlp_block !== nothing
        num_con += length(nlp_block.constraint_bounds)
    end
    g_l = fill(-Inf, num_con)
    g_u = fill(Inf, num_con)
    if nlp_block !== nothing
        for (i, bound) in enumerate(nlp_block.constraint_bounds)
            g_l[i] = bound.lower
            g_u[i] = bound.upper
        end
    end

    # ==========================================================================
    # Intialize the variables. We need to form a mapping between the MOI
    # VariableIndex and an Int in order to replace instances of
    # `x[VariableIndex]` with `x[i]` in the expression graphs.
    variables = MOI.get(model, MOI.ListOfVariableIndices())
    num_var = length(variables)
    variable_map = Dict{MOI.VariableIndex, Int}()
    for (i, variable) in enumerate(variables)
        variable_map[variable] = i
    end

    # ==========================================================================
    # Extract variable bounds.
    x_l = fill(-Inf, num_var)
    x_u = fill(Inf, num_var)
    for set_type in MOI_SCALAR_SETS
        for c_ref in MOI.get(model,
            MOI.ListOfConstraintIndices{MOI.SingleVariable, set_type}())
            c_func = MOI.get(model, MOI.ConstraintFunction(), c_ref)
            c_set = MOI.get(model, MOI.ConstraintSet(), c_ref)
            v_index = variable_map[c_func.variable]
            lower, upper = set_to_bounds(c_set)
            # Note that there might be multiple bounds on the same variable
            # (e.g., LessThan and GreaterThan), so we should only update bounds
            # if they differ from the defaults.
            if lower > -Inf
                x_l[v_index] = lower
            end
            if upper < Inf
                x_u[v_index] = upper
            end
        end
    end

    # ==========================================================================
    # We have to convert all ScalarAffineFunction-in-Set constraints to an
    # expression graph.
    scalar_constraint_expr = Expr[]
    for func_type in MOI_SCALAR_FUNCTIONS
        for set_type in MOI_SCALAR_SETS
            for c_ref in MOI.get(model, MOI.ListOfConstraintIndices{
                    func_type, set_type}())
                c_func = MOI.get(model, MOI.ConstraintFunction(), c_ref)
                c_set = MOI.get(model, MOI.ConstraintSet(), c_ref)
                expr = moi_to_expr_graph(c_func, c_set, variable_map)
                push!(scalar_constraint_expr, expr)
                lower, upper = set_to_bounds(c_set)
                push!(g_l, lower)
                push!(g_u, upper)
            end
        end
    end

    # ==========================================================================
    # MOI objective
    obj_type = MOI.get(model, MOI.ObjectiveFunctionType())
    obj_func = MOI.get(model, MOI.ObjectiveFunction{obj_type}())
    obj_func_expr = func_to_expr_graph(obj_func, variable_map)
    if obj_func_expr == :(+ 0.0)
        obj_func_expr = nothing
    end

    # ==========================================================================
    # Build the nlp_evaluator
    nlp_evaluator = NLPEvaluator(moi_nlp_evaluator, variable_map, num_con,
        obj_func_expr, scalar_constraint_expr)

    # ==========================================================================
    # Create the MathProgBase model. Note that we pass `num_con` and the number
    # of linear constraints.
    mpb_model = MPB.NonlinearModel(mpb_solver)
    MPB.loadproblem!(mpb_model, num_var,
        num_con + length(scalar_constraint_expr), x_l, x_u, g_l, g_u, sense,
        nlp_evaluator)

    # ==========================================================================
    # Set any variables to :Bin if they are in ZeroOne and :Int if they are
    # Integer. The default is just :Cont.
    x_cat = fill(:Cont, num_var)
    for set_type in (MOI.ZeroOne, MOI.Integer)
        for c_ref in MOI.get(model,
            MOI.ListOfConstraintIndices{MOI.SingleVariable, set_type}())
            c_func = MOI.get(model, MOI.ConstraintFunction(), c_ref)
            c_set = MOI.get(model, MOI.ConstraintSet(), c_ref)
            v_index = variable_map[c_func.variable]
            x_cat[v_index] = set_to_cat(c_set)
        end
    end
    MPB.setvartype!(mpb_model, x_cat)

    # ==========================================================================
    # Set the VariablePrimalStart attributes for variables.
    variable_primal_start = fill(0.0, num_var)
    for (i, variable) in enumerate(variables)
        try
            start_val = MOI.get(model, MOI.VariablePrimalStart(), variable)
            if start_val !== nothing
                variable_primal_start[i] = start_val
            end
        catch
        end
    end
    MPB.setwarmstart!(mpb_model, variable_primal_start)

    # ==========================================================================
    # Set the VariablePrimalStart attributes for variables.
    MPB.optimize!(mpb_model)

    # ==========================================================================
    # Extract and save the MathProgBase solution.
    primal_solution = Dict{MOI.VariableIndex, Float64}()
    for (variable, sol) in zip(variables, MPB.getsolution(mpb_model))
        primal_solution[variable] = sol
    end
    mpb_solution = MPBSolution(
        MPB.status(mpb_model),
        MPB.getobjval(mpb_model),
        primal_solution
    )
    MOI.set(model, MPBSolutionAttribute(), mpb_solution)
    return
end

# ==============================================================================
# MOI accessors for the solution info.

function MOI.get(model::Model, ::MOI.VariablePrimal, var::MOI.VariableIndex)
    mpb_solution = try
        MOI.get(model, MPBSolutionAttribute())
    catch
        nothing
    end
    if mpb_solution === nothing
        return nothing
    end
    return mpb_solution.primal_solution[var]
end

function MOI.get(model::Model, ::MOI.ConstraintPrimal, idx::MOI.ConstraintIndex)
    return MOIU.get_fallback(model, MOI.ConstraintPrimal(), idx)
end

function MOI.get(model::Model, ::MOI.ObjectiveValue)
    mpb_solution = try
        MOI.get(model, MPBSolutionAttribute())
    catch
        nothing
    end
    if mpb_solution === nothing
        return nothing
    end
    return mpb_solution.objective_value
end

function MOI.get(model::Model, ::MOI.TerminationStatus)
    mpb_solution = try
        MOI.get(model, MPBSolutionAttribute())
    catch
        nothing
    end
    if mpb_solution === nothing
        return MOI.OPTIMIZE_NOT_CALLED
    end
    status = mpb_solution.status
    if status == :Optimal
        # TODO(odow): this is not always the case. What if Ipopt solves a
        # convex problem?
        return MOI.LOCALLY_SOLVED
    elseif status == :Infeasible
        return MOI.INFEASIBLE
    elseif status == :Unbounded
        return MOI.DUAL_INFEASIBLE
    elseif status == :UserLimit
        return MOI.OTHER_LIMIT
    elseif status == :Error
        return MOI.OTHER_ERROR
    end
    return MOI.OTHER_ERROR
end

function MOI.get(model::Model, ::MOI.PrimalStatus)
    mpb_solution = try
        MOI.get(model, MPBSolutionAttribute())
    catch
        nothing
    end
    if mpb_solution === nothing
        return MOI.NO_SOLUTION
    end
    status = mpb_solution.status
    if status == :Optimal
        return MOI.FEASIBLE_POINT
    end
    return MOI.NO_SOLUTION
end

function MOI.get(model::Model, ::MOI.DualStatus)
    return MOI.NO_SOLUTION
end

function MOI.get(model::Model, ::MOI.ResultCount)
    if MOI.get(model, MOI.PrimalStatus()) == MOI.FEASIBLE_POINT
        return 1
    else
        return 0
    end
end

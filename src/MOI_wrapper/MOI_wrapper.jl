"""
MOI_wrapper.jl defines the Juniper.Optimizer struct
with all mandatory MOI functions overloaded
"""

"""
Optimizer struct
"""
mutable struct Optimizer <: MOI.AbstractOptimizer
    inner::Union{JuniperProblem, Nothing}
    model_cache::MOIU.UniversalFallback{MOIU.Model{Float64}}
    options::SolverOptions
end

MOI.is_valid(model::Optimizer, index::MOI.Index) = MOI.is_valid(model.model_cache, index)

MOI.get(::Optimizer, ::MOI.SolverName) = "Juniper"
MOI.get(::Optimizer, ::MOI.SolverVersion) = "v0.7.0"

MOI.supports(::Optimizer, ::MOI.Silent) = true
MOI.supports(::Optimizer, ::MOI.NumberOfThreads) = true
MOI.supports(::Optimizer, ::MOI.TimeLimitSec) = true
MOI.supports(::Optimizer, ::MOI.RawOptimizerAttribute) = true

function MOI.set(model::Optimizer, ::MOI.Silent, value::Bool)
    if value
        model.options.log_levels = []
    end
    model.options.silent = value
    return
end

function MOI.set(model::Optimizer, ::MOI.NumberOfThreads, value::Union{Nothing,Int})
    if value === nothing
        model.options.processors = 1
    else
        model.options.processors = value
    end
    return
end

function MOI.set(model::Optimizer, ::MOI.TimeLimitSec, value::Union{Nothing,Float64})
    if value === nothing
        model.options.time_limit = Inf
    else
        model.options.time_limit = value
    end
    return
end

function MOI.set(model::Optimizer, p::MOI.RawOptimizerAttribute, value)
    p_symbol = Symbol(p.name)
    if in(p_symbol, fieldnames(SolverOptions))
        type_of_param = fieldtype(SolverOptions, p_symbol)
        if hasmethod(convert, (Type{type_of_param}, typeof(value)))
            passed_checks = true
            if p_symbol == :traverse_strategy
                if !(value in [:BFS, :DFS, :DBFS])
                    passed_checks = false
                    @error "Traverse strategy $(value) is not supported. Use one of `[:BFS, :DFS, :DBFS]`."
                end
            end
            if p_symbol == :branch_strategy
                if !(value in [:StrongPseudoCost, :PseudoCost, :Reliability, :MostInfeasible])
                    passed_checks = false
                    @error "Branch strategy $(value) is not supported. Use one of `[:StrongPseudoCost, :PseudoCost, :Reliability, :MostInfeasible]`."
                end
            end
            passed_checks && setfield!(model.options, p_symbol, convert(type_of_param, value))
        else
            @error "The option $(p.name) has a different type ($(type_of_param))"
        end
    else
        @error "The option $(p.name) doesn't exist."
    end
    return
end

# Returns the number of solutions that can be retrieved
MOI.get(model::Optimizer, ::MOI.ResultCount) = length(model.inner.solutions)

MOI.get(model::Optimizer, ::MOI.NumberOfThreads) = model.options.processors

MOI.get(model::Optimizer, ::MOI.TimeLimitSec) = model.options.time_limit

MOI.get(model::Optimizer, ::MOI.Silent) = model.options.silent

function MOI.get(model::Optimizer, p::MOI.RawOptimizerAttribute)
    p_symbol = Symbol(p.name)
    if in(p_symbol, fieldnames(SolverOptions))
        return getfield(model.options, p_symbol)
    end
    @error "The option $(p.name) doesn't exist."
end

function register(s::Symbol, dimension::Integer, f::Function; autodiff::Bool=false)
    return RegisteredFunction(s, dimension, f, nothing, nothing, autodiff)
end

function register(s::Symbol, dimension::Integer, f::Function, gradf::Function; autodiff::Bool=false)
    return RegisteredFunction(s, dimension, f, gradf, nothing, autodiff)
end

function register(s::Symbol, dimension::Integer, f::Function, gradf::Function, grad2f::Function; autodiff::Bool=false)
    return RegisteredFunction(s, dimension, f, gradf, grad2f, autodiff)
end

"""
Optimizer struct constructor
"""
function Optimizer(;options...)

    solver_options = combine_options(options)

    return Optimizer(
    nothing,
    MOIU.UniversalFallback(MOIU.Model{Float64}()),
    solver_options)
end

function Optimizer(options::Vector{Pair{String,Any}})
    symbol_options = Dict{Symbol, Any}()
    for option in options
        symbol_options[Symbol(option.first)] = option.second
    end
    return Optimizer(symbol_options)
end

function Optimizer(options::Dict{Symbol,Any})

    solver_options = combine_options(options)

    return Optimizer(
    nothing,
    MOIU.UniversalFallback(MOIU.Model{Float64}()),
    solver_options)
end

"""
Printing the optimizer
"""
function Base.show(io::IO, model::Optimizer)
    println("A Juniper MathOptInterface model with backend")
    return
end

"""
Copy constructor for the optimizer
"""
MOIU.supports_default_copy_to(model::Optimizer, copy_names::Bool) = true
function MOI.copy_to(model::Optimizer, src::MOI.ModelLike; kws...)
    return MOI.copy_to(model.model_cache, src; kws...)
end

"""
``MOI.is_empty(model::Optimizer)`` overload for Alpine.Optimizer
"""
MOI.is_empty(model::Optimizer) = MOI.is_empty(model.model_cache)

"""
``MOI.empty!(model::Optimizer)`` overload for Alpine.Optimizer
"""
function MOI.empty!(model::Optimizer)
    model.inner = nothing
    MOI.empty!(model.model_cache)
end

"""
ordering of constraints provided to Juniper.jl
"""
linear_le_offset(model::Optimizer) = 0
linear_ge_offset(model::Optimizer) = length(model.linear_le_constraints)
linear_eq_offset(model::Optimizer) = linear_ge_offset(model) + length(model.linear_ge_constraints)
quadratic_le_offset(model::Optimizer) = linear_eq_offset(model) + length(model.linear_eq_constraints)
quadratic_ge_offset(model::Optimizer) = quadratic_le_offset(model) + length(model.quadratic_le_constraints)
quadratic_eq_offset(model::Optimizer) = quadratic_ge_offset(model) + length(model.quadratic_ge_constraints)
nlp_constraint_offset(model::Optimizer) = quadratic_eq_offset(model) + length(model.quadratic_eq_constraints)

function replace_solution!(m::JuniperProblem, tree::BnBTreeObj)
    status_dict = Dict{Symbol,MOI.TerminationStatusCode}()
    status_dict[:Time] = MOI.TIME_LIMIT
    status_dict[:Infeasible] = !tree.global_solver ? MOI.LOCALLY_INFEASIBLE : MOI.INFEASIBLE
    status_dict[:MipGap] = MOI.OBJECTIVE_LIMIT
    status_dict[:BestObjStop] = MOI.OBJECTIVE_LIMIT
    status_dict[:EnoughSolutions] = MOI.SOLUTION_LIMIT

    if isdefined(tree, :incumbent)
        incumbent = tree.incumbent
        m.objval = incumbent.objval
        m.solution = incumbent.solution
    end
    if tree.limit == :None
        if tree.incumbent.only_almost
            if !tree.global_solver
                m.status = MOI.ALMOST_LOCALLY_SOLVED
            else
                m.status = MOI.ALMOST_OPTIMAL
            end
        else
            if !tree.global_solver
                m.status = MOI.LOCALLY_SOLVED
            else
                m.status = MOI.OPTIMAL
            end
        end
    else
        m.status = status_dict[tree.limit]
    end
    m.best_bound = tree.best_bound
end

"""
``MOI.optimize!()`` for Juniper
"""
function MOI.optimize!(model::Optimizer)
    # feasibility_pump only if mip solver exists
    if model.options.feasibility_pump && model.options.mip_solver === nothing
        model.options.feasibility_pump = false
    end
    Random.seed!(JUNIPER_RNG, model.options.seed)

    # fill JuniperProblem
    model.inner = JuniperProblem()
    @views jp = model.inner
    init_juniper_problem!(jp, model)

    ps = jp.options.log_levels
    jp.debugDict = Dict{Any,Any}()

    (:All in ps || :AllOptions in ps) && print_options(jp;all=true)
    (:Options in ps) && print_options(jp;all=false)

    if !jp.options.fixed_gain_mu && jp.obj_sense == :Max
        jp.options.gain_mu = 1-jp.options.gain_mu
    end

    nw = nworkers()
    if nw < jp.options.processors
        jp.options.processors = nw
        @warn "Julia was started with less processors than you defined in your options. Start julia with: `julia -p "*string(jp.options.processors)*"`"
    end
    # set incumbent to nothing might be updated using start values or the feasibility_pump
    incumbent = nothing

    create_root_model!(model, jp)
    (:All in ps || :Info in ps) && print_info(jp)
    # fix primal start to check if we have an incumbent
    fix_primal_start!(jp)

    incumbent = solve_root_incumbent_model(jp)

    unfix_primal_start!(jp)
    relax_start_time = time()
    restarts = solve_root_model!(jp)
    jp.relaxation_time = time()-relax_start_time

    (:All in ps || :Info in ps) && println("Status of relaxation: ", jp.relaxation_status)
    jp.soltime = time()-jp.start_time

    jp.options.debug && debug_fill_basic(jp.debugDict,jp,restarts)

    # if infeasible or unbounded => return
    if !state_is_optimal(jp.relaxation_status; allow_almost=jp.options.allow_almost_solved)
        jp.status = jp.relaxation_status
        if jp.options.debug && jp.options.debug_write
            write(jp.options.debug_file_path, JSON.json(jp.debugDict))
        end
        return
    end

    (:All in ps || :Info in ps || :Timing in ps) && println("Time for relaxation: ", jp.soltime)

    jp.relaxation_objval   = MOI.get(jp.model, MOI.ObjectiveValue())
    jp.relaxation_solution = MOI.get(jp.model, MOI.VariablePrimal(), jp.x)

    jp.options.debug && debug_objective(jp.debugDict,jp)
    # TODO free model for Knitro

    (:All in ps || :Info in ps || :Timing in ps) && println("Relaxation Obj: ", jp.relaxation_objval)


    only_almost_solved = false
    if jp.num_disc_var > 0
        if jp.options.feasibility_pump
            fpump_incumbent = fpump(model,jp)
            if fpump_incumbent !== nothing
                if incumbent === nothing
                    incumbent = fpump_incumbent
                else
                    factor = jp.obj_sense == :Min ? -1 : 1
                    # if found better incumbent
                    if factor*fpump_incumbent.objval > factor*incumbent.objval
                        incumbent = fpump_incumbent
                    end
                end
            end
        end
        bnbtree = init(jp.start_time, jp; incumbent = incumbent)
        solvemip(bnbtree)

        replace_solution!(jp, bnbtree)
        jp.nsolutions = bnbtree.nsolutions
    else
        jp.nsolutions = 1
        jp.best_bound = try
            MOI.get(jp.model, MOI.ObjectiveBound())
        catch
            MOI.get(jp.model, MOI.ObjectiveValue())
        end
        jp.status = jp.relaxation_status
        jp.objval = jp.relaxation_objval
        jp.solution = jp.relaxation_solution
    end
    jp.soltime = time()-jp.start_time

    (:All in ps || :Info in ps) && println("Obj: ",jp.objval)

    if length(jp.solutions) == 0
        push!(jp.solutions, SolutionObj(jp.solution, jp.objval))
    end

    jp.options.debug && debug_set_solution(jp.debugDict,jp)
    if jp.options.debug && jp.options.debug_write
        write(jp.options.debug_file_path, JSON.json(jp.debugDict))
    end

end

getnsolutions(m::JuniperProblem) = m.nsolutions
getsolutions(m::JuniperProblem) = m.solutions
getnbranches(m::JuniperProblem) = m.nbranches

# `UniversalFallback` supports everything so we can return `true`.
MOI.supports(model::Optimizer, attr::MOI.AbstractModelAttribute) = true
MOI.set(model::Optimizer, attr::MOI.AbstractModelAttribute, value) = MOI.set(model.model_cache, attr, value)
MOI.get(model::Optimizer, attr::MOI.AbstractModelAttribute) = MOI.get(model.model_cache, attr)

"""
MOI variables
"""

MOI.add_variable(model::Optimizer) = MOI.add_variable(model.model_cache)
MOI.add_variables(model::Optimizer, n) = MOI.add_variables(model.model_cache, n)

# `UniversalFallback` supports everything so we can return `true`.
MOI.supports(::Optimizer, ::MOI.AbstractVariableAttribute, ::Type{MOI.VariableIndex}) = true
function MOI.set(model::Optimizer, attr::MOI.AbstractVariableAttribute, vi::MOI.VariableIndex, value)
    MOI.set(model.model_cache, attr, vi, value)
    return
end
function MOI.get(model::Optimizer, attr::MOI.AbstractVariableAttribute, vi::MOI.VariableIndex)
    return MOI.get(model.model_cache, attr, vi)
end

"""
MOI constraints
"""

# `UniversalFallback` supports everything so we can return `true`.
MOI.supports_constraint(::Optimizer, ::Type{<:MOI.AbstractFunction}, ::Type{<:MOI.AbstractSet}) = true
function MOI.add_constraint(model::Optimizer, func::MOI.AbstractFunction, set::MOI.AbstractSet)
    return MOI.add_constraint(model.model_cache, func, set)
end

# `UniversalFallback` supports everything so we can return `true`.
MOI.supports(::Optimizer, ::MOI.AbstractConstraintAttribute, ::Type{<:MOI.ConstraintIndex}) = true
function MOI.set(model::Optimizer, attr::MOI.AbstractConstraintAttribute, ci::MOI.ConstraintIndex, value)
    MOI.set(model.model_cache, attr, ci, value)
    return
end
function MOI.get(model::Optimizer, attr::MOI.AbstractConstraintAttribute, ci::MOI.ConstraintIndex)
    return MOI.get(model.model_cache, attr, ci)
end

function MOI.get(
    model::Optimizer,
    attr::Union{MOI.ConstraintFunction,MOI.ConstraintSet},
    ci::MOI.ConstraintIndex,
)
    return MOI.get(model.model_cache, attr, ci)
end

include("results.jl")

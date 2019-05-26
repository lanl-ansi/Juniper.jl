"""
MOI_wrapper.jl defines the Juniper.Optimizer struct
with all mandatory MOI functions overloaded
"""


""" 
MOI functions, sets and, other type definitions
"""
# indices
const VI = MOI.VariableIndex
const CI = MOI.ConstraintIndex

# sets
const BOUNDS = Union{
    MOI.EqualTo{Float64}, 
    MOI.GreaterThan{Float64},
    MOI.LessThan{Float64}, 
    MOI.Interval{Float64}
}

const VAR_TYPES = Union{
    MOI.ZeroOne, 
    MOI.Integer
}

# other MOI types
const AFF_TERM = MOI.ScalarAffineTerm{Float64}
const QUAD_TERM = MOI.ScalarQuadraticTerm{Float64}

"""
Variable information struct definition 
""" 
mutable struct VariableInfo
    lower_bound::Float64  # May be -Inf even if has_lower_bound == true
    has_lower_bound::Bool # false implies lower_bound == -Inf
    upper_bound::Float64  # May be Inf even if has_upper_bound == true
    has_upper_bound::Bool # false implies upper_bound == Inf
    is_fixed::Bool        # Implies lower_bound == upper_bound and !has_lower_bound and !has_upper_bound
    is_binary::Bool       # Implies lower_bound == 0, upper_bound == 1 and is MOI.ZeroOne
    is_integer::Bool      # Implies variable is MOI.Integer
    start::Real           # Primal start
    name::String
end
VariableInfo() = VariableInfo(-Inf, false, Inf, false, false, false, false, 0.0, "")

"""
Optimizer struct
"""  
mutable struct Optimizer <: MOI.AbstractOptimizer
    inner::Union{JuniperProblem, Nothing}
    variable_info::Vector{VariableInfo}
    nlp_data::MOI.NLPBlockData
    sense::MOI.OptimizationSense 
    objective::Union{SVF, SAF, SQF, Nothing}
    linear_le_constraints::Vector{Tuple{SAF, MOI.LessThan{Float64}}}
    linear_ge_constraints::Vector{Tuple{SAF, MOI.GreaterThan{Float64}}}
    linear_eq_constraints::Vector{Tuple{SAF, MOI.EqualTo{Float64}}}
    quadratic_le_constraints::Vector{Tuple{SQF, MOI.LessThan{Float64}}}
    quadratic_ge_constraints::Vector{Tuple{SQF, MOI.GreaterThan{Float64}}}
    quadratic_eq_constraints::Vector{Tuple{SQF, MOI.EqualTo{Float64}}}
    options::SolverOptions
end

MOI.get(::Optimizer, ::MOI.SolverName) = "Juniper"

"""
EmptyNLPEvaluator struct and associated functions 
"""
struct EmptyNLPEvaluator <: MOI.AbstractNLPEvaluator end
MOI.features_available(::EmptyNLPEvaluator) = [:Grad, :Jac, :Hess]
MOI.initialize(::EmptyNLPEvaluator, features) = nothing
MOI.eval_objective(::EmptyNLPEvaluator, x) = NaN
function MOI.eval_constraint(::EmptyNLPEvaluator, g, x)
    @assert length(g) == 0
    return
end
function MOI.eval_objective_gradient(::EmptyNLPEvaluator, g, x)
    fill!(g, 0.0)
    return
end
MOI.jacobian_structure(::EmptyNLPEvaluator) = Tuple{Int64,Int64}[]
MOI.hessian_lagrangian_structure(::EmptyNLPEvaluator) = Tuple{Int64,Int64}[]
function MOI.eval_constraint_jacobian(::EmptyNLPEvaluator, J, x)
    @assert length(J) == 0
    return
end
function MOI.eval_hessian_lagrangian(::EmptyNLPEvaluator, H, x, σ, μ)
    @assert length(H) == 0
    return
end
empty_nlp_data() = MOI.NLPBlockData([], EmptyNLPEvaluator(), false)

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
    [], 
    empty_nlp_data(), 
    MOI.FEASIBILITY_SENSE, 
    nothing, 
    [], [], [], # linear constraints 
    [], [], [], # quadratic constraints
    solver_options)
end 

function Optimizer(options::Dict{Symbol,Any}) 
    
    solver_options = combine_options(options)

    return Optimizer(
    nothing, 
    [], 
    empty_nlp_data(), 
    MOI.FEASIBILITY_SENSE, 
    nothing, 
    [], [], [], # linear constraints 
    [], [], [], # quadratic constraints
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
    return MOI.Utilities.automatic_copy_to(model, src; kws...)
end

"""
``MOI.is_empty(model::Optimizer)`` overload for Alpine.Optimizer
"""
function MOI.is_empty(model::Optimizer)
    return isempty(model.variable_info) &&
        model.nlp_data.evaluator isa EmptyNLPEvaluator &&
        model.sense == MOI.FEASIBILITY_SENSE &&
        isempty(model.linear_le_constraints) &&
        isempty(model.linear_ge_constraints) &&
        isempty(model.linear_eq_constraints) &&
        isempty(model.quadratic_le_constraints) &&
        isempty(model.quadratic_ge_constraints) &&
        isempty(model.quadratic_eq_constraints)
end

"""
``MOI.empty!(model::Optimizer)`` overload for Alpine.Optimizer
"""
function MOI.empty!(model::Optimizer)
    model.inner = nothing
    empty!(model.variable_info)
    model.nlp_data = empty_nlp_data()
    model.sense = MOI.FEASIBILITY_SENSE
    model.objective = nothing
    empty!(model.linear_le_constraints)
    empty!(model.linear_ge_constraints)
    empty!(model.linear_eq_constraints)
    empty!(model.quadratic_le_constraints)
    empty!(model.quadratic_ge_constraints)
    empty!(model.quadratic_eq_constraints)
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


function info_array_of_variables(variable_info::Vector{VariableInfo}, attr::Symbol)
    len_var_info = length(variable_info)
    type_dict = get_type_dict(variable_info[1])
    result = Array{type_dict[attr], 1}(undef, len_var_info)
    for i = 1:len_var_info
        result[i] = getfield(variable_info[i], attr)
        # if type is binary then set bounds correctly
        if result[i] < 0 && attr == :lower_bound && getfield(variable_info[i], :is_binary)
            result[i] = 0
        end
        if result[i] > 1 && attr == :upper_bound && getfield(variable_info[i], :is_binary)
            result[i] = 1
        end
    end
    return result
end

function replace_solution!(m::JuniperProblem, best_known)
    m.solution = best_known.solution
    m.objval = best_known.objval
    m.status = best_known.status
    m.best_bound = best_known.best_bound # is reasonable for gap or time limit
end

"""
``MOI.optimize!()`` for Juniper
""" 
function MOI.optimize!(model::Optimizer)
    Random.seed!(1)
    MOI.initialize(model.nlp_data.evaluator, [:ExprGraph])
    
    if ~isa(model.nlp_data.evaluator, EmptyNLPEvaluator)

    else 
        @info "no explicit NLP constraints or objective provided using @NLconstraint or @NLobjective macros"
    end 

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

    create_root_model!(model, jp)
    relax_start_time = time()
    restarts = solve_root_model!(jp)
    jp.relaxation_time = time()-relax_start_time

    (:All in ps || :Info in ps) && println("Status of relaxation: ", jp.status)
    jp.soltime = time()-jp.start_time

    jp.options.debug && debug_fill_basic(jp.debugDict,jp,restarts)

    # if infeasible or unbounded => return
    if !state_is_optimal(jp.status; allow_almost=true)
        if jp.options.debug && jp.options.debug_write
            write(jp.options.debug_file_path, JSON.json(jp.debugDict))
        end
        return jp.status
    end

    (:All in ps || :Info in ps || :Timing in ps) && println("Time for relaxation: ", jp.soltime)
    
    backend     = JuMP.backend(jp.model)
    jp.objval   = JuMP.objective_value(jp.model)
    jp.solution = JuMP.value.(jp.x)

    jp.options.debug && debug_objective(jp.debugDict,jp)
    # TODO free model for Knitro

    (:All in ps || :Info in ps || :Timing in ps) && println("Relaxation Obj: ", jp.objval)

    # set incumbent to nothing might be updated using the feasibility_pump
    inc_sol, inc_obj = nothing, nothing
    only_almost_solved = false
    if jp.num_disc_var > 0
        if jp.options.feasibility_pump
            inc_sol, inc_obj, only_almost_solved = fpump(model,jp)
        end
        bnbtree = init(jp.start_time, jp; inc_sol = inc_sol, inc_obj = inc_obj, only_almost_solved = only_almost_solved)
        best_known = solvemip(bnbtree)

        replace_solution!(jp, best_known)
        jp.nsolutions = bnbtree.nsolutions
    else
        jp.nsolutions = 1
        jp.best_bound = try 
            JuMP.objective_bound(jp.model)
        catch
            JuMP.objective_value(jp.model)
        end
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

include("variables.jl")
include("constraints.jl")
include("objective.jl")
include("results.jl")
include("nlp.jl")
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

# functions
const SVF = MOI.SingleVariable
const SAF = MOI.ScalarAffineFunction{Float64}
const SQF = MOI.ScalarQuadraticFunction{Float64}
const VECTOR = MOI.VectorOfVariables

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
    name::String
end
VariableInfo() = VariableInfo(-Inf, false, Inf, false, false, false, false, "")

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

"""
Printing the optimizer 
"""
function Base.show(io::IO, model::Optimizer)
    println(io, "A MathOptInterface model with backend:")
    println(io, model.inner)
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
    println("Type of attr (",attr,") is: ", type_dict[attr])
    result = Array{type_dict[attr], 1}(undef, len_var_info)
    for i = 1:len_var_info
        result[i] = getfield(variable_info[i], attr)
    end
    return result
end

"""
``MOI.optimize!()`` for Juniper
""" 
function MOI.optimize!(model::Optimizer)   
    num_variables = length(model.variable_info)
    num_linear_le_constraints = length(model.linear_le_constraints)
    num_linear_ge_constraints = length(model.linear_ge_constraints)
    num_linear_eq_constraints = length(model.linear_eq_constraints)
    num_quadratic_le_constraints = length(model.quadratic_le_constraints)
    num_quadratic_ge_constraints = length(model.quadratic_ge_constraints)
    num_quadratic_eq_constraints = length(model.quadratic_eq_constraints)

    println("linear eq constraints: ")
    println(model.linear_eq_constraints)
    println("Objective sense: ")
    println(model.sense)
    println("Objective: ")
    println(model.objective)

    
    if ~isa(model.nlp_data.evaluator, EmptyNLPEvaluator)

    else 
        @info "no explicit NLP constraints or objective provided using @NLconstraint or @NLobjective macros"
    end 

    # fill JuniperProblem
    model.inner = JuniperProblem()
    @views jp = model.inner
    jp.start_time = time()

    jp.nl_solver = model.options.nl_solver
    if model.options.mip_solver != nothing
        jp.mip_solver = model.options.mip_solver
    end
    jp.options = model.options    
    if model.sense == MOI.MIN_SENSE 
        jp.obj_sense = :Min
    else
        jp.obj_sense = :Max
    end
    jp.l_var = info_array_of_variables(model.variable_info, :lower_bound)
    jp.u_var = info_array_of_variables(model.variable_info, :upper_bound)
    jp.num_var = length(model.variable_info)

    
    ps = jp.options.log_levels
    println("ps: ", ps)
    jp.debugDict = Dict{Any,Any}()
    if !jp.options.fixed_gain_mu && jp.obj_sense == :Max
        jp.options.gain_mu = 1-jp.options.gain_mu
    end
    (:All in ps || :AllOptions in ps) && print_options(jp;all=true)
    (:Options in ps) && print_options(jp;all=false)

    nw = nworkers()
    if nw < jp.options.processors
        jp.options.processors = nw
        @warn "Julia was started with less processors then you define in your options"
    end

    create_root_model!(model, jp)
    relax_start_time = time()
    restarts = solve_root_model!(jp)
    jp.relaxation_time = time()-relax_start_time

    (:All in ps || :Info in ps) && println("Status of relaxation: ", jp.status)
    jp.soltime = time()-jp.start_time

    jp.options.debug && debug_fill_basic(jp.debugDict,jp,restarts)

    # if infeasible or unbounded => return
    if jp.status != MOI.OPTIMAL && jp.status != MOI.LOCALLY_SOLVED
        if jp.options.debug && jp.options.debug_write
            write(jp.options.debug_file_path, JSON.json(jp.debugDict))
        end
        return jp.status
    end

    (:All in ps || :Info in ps || :Timing in ps) && println("Time for relaxation: ", jp.soltime)
    
    backend     = JuMP.backend(jp.model)
    # TODO: doesn't work atm
    jp.objval   = JuMP.objective_value(jp.model)
    jp.solution = get_primal_values(backend)
    println("Objval: ", jp.objval)
    println("solution: ", jp.solution)

    jp.options.debug && debug_objective(jp.debugDict,m)
    

end 

include("variables.jl")
include("constraints.jl")
include("objective.jl")
include("results.jl")
include("nlp.jl")
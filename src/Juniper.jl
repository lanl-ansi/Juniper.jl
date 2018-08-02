module Juniper

using MathProgBase
using JuMP

# Options for the solver (more details like defaults in solver.jl)
type SolverOptions
    log_levels                          :: Vector{Symbol}
    atol                                :: Float64
    num_resolve_root_relaxation         :: Int64
    branch_strategy                     :: Symbol
    gain_mu                             :: Float64
    strong_branching_perc               :: Float64
    strong_branching_nsteps             :: Int64
    strong_branching_approx_time_limit  :: Float64
    strong_restart                      :: Bool
    reliability_branching_threshold     :: Int64
    reliability_branching_perc          :: Float64
    incumbent_constr                    :: Bool
    obj_epsilon                         :: Float64
    time_limit                          :: Float64
    mip_gap                             :: Float64
    best_obj_stop                       :: Float64
    solution_limit                      :: Int64
    all_solutions                       :: Bool
    list_of_solutions                   :: Bool
    processors                          :: Int64
    traverse_strategy                   :: Symbol
    feasibility_pump                    :: Bool
    feasibility_pump_time_limit         :: Float64
    feasibility_pump_tolerance_counter  :: Int64
    tabu_list_length                    :: Int64
    num_resolve_nlp_feasibility_pump    :: Int64
    mip_solver                          :: Union{Void, MathProgBase.AbstractMathProgSolver}
    
    # only for testing
    force_parallel                      :: Bool

    fixed_gain_mu                       :: Bool
end

function Base.show(io::IO, opts::SolverOptions) 
    longest_field_name = maximum([length(string(fname)) for fname in fieldnames(SolverOptions)])+2
    for name in fieldnames(SolverOptions)
        sname = string(name)
        pname = sname*repeat(" ", longest_field_name-length(sname))
        println(io, pname, ": ", getfield(opts,name))
    end
end

include("solver.jl")
include("util.jl")
include("model.jl")
include("BnBTree.jl")


end

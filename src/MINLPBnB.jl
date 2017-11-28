module MINLPBnB

using MathProgBase
using JuMP

# Options for the solver (more details like defaults in solver.jl)
type SolverOptions
    log_levels                  :: Vector{Symbol}
    branch_strategy             :: Symbol
    gain_mu                     :: Float64
    strong_branching_perc       :: Float64
    strong_branching_nsteps     :: Int64
    strong_restart              :: Bool
    incumbent_constr            :: Bool
    obj_epsilon                 :: Float64
    time_limit                  :: Float64
    mip_gap                     :: Float64
    best_obj_stop               :: Float64
    solution_limit              :: Int64
    all_solutions               :: Bool
    list_of_solutions           :: Bool
    processors                  :: Int64
    traverse_strategy           :: Symbol
    feasibility_pump            :: Bool
    feasibility_pump_time_limit :: Int64
    mip_solver                  :: Union{Void,MathProgBase.AbstractMathProgSolver} # TODO: Maybe not union?
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
include("model.jl")
include("BnBTree.jl")

end

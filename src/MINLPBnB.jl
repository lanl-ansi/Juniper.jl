module MINLPBnB

using MathProgBase
using JuMP

# Options for the solver (more details like defaults in solver.jl)
type SolverOptions
    log_levels                  :: Vector{Symbol}
    branch_strategy             :: Symbol
    strong_branching_nvars      :: Int64
    strong_branching_nsteps     :: Int64
    strong_restart              :: Bool
    incumbent_constr            :: Bool
    time_limit                  :: Float64
    mip_gap                     :: Float64
    best_obj_stop               :: Float64
    solution_limit              :: Int64
    all_solutions               :: Bool
    list_of_solutions           :: Bool
end

include(Pkg.dir("MINLPBnB")*"/src/solver.jl")
include(Pkg.dir("MINLPBnB")*"/src/model.jl")
include(Pkg.dir("MINLPBnB")*"/src/BnBTree.jl")

end

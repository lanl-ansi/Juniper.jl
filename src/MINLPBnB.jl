module MINLPBnB

using MathProgBase
using JuMP

type SolverOptions
    log_levels                  :: Vector{Symbol}
    branch_strategy             :: Symbol
    strong_branching_nvars      :: Int64
    strong_branching_nsteps     :: Int64
    strong_restart              :: Bool
    incumbent_constr            :: Bool
    time_limit                  :: Float64
    mip_gap                     :: Float64
end

include("solver.jl")
include("model.jl")
include("BnBTree.jl")

end

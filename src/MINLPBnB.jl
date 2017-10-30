module MINLPBnB

using MathProgBase
using JuMP

type SolverOptions
    log_levels                  :: Vector{Symbol}
    branch_strategy             :: Symbol
    strong_branching_nvars      :: Int64
    strong_branching_nlevels    :: Int64
    strong_restart              :: Bool
end

include("solver.jl")
include("model.jl")
include("BnBTree.jl")

end

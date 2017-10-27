module MINLPBnB

using MathProgBase
using JuMP

mutable struct SolverOptions
    log_levels                  :: Vector{Symbol}
    branch_strategy             :: Symbol
    strong_branching_nvars      :: Int64
    strong_branching_nlevels    :: Int64
end

include("solver.jl")
include("model.jl")
include("BnBTree.jl")

end

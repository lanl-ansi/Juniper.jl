module Juniper

using JSON
using LinearAlgebra
using Random
const JUNIPER_RNG = MersenneTwister(1)

using Distributed
using Statistics

import MutableArithmetics
const MA = MutableArithmetics

using MathOptInterface
const MOI = MathOptInterface
const MOIU = MOI.Utilities

# functions
const SAF = MOI.ScalarAffineFunction{Float64}
const SQF = MOI.ScalarQuadraticFunction{Float64}
const VECTOR = MOI.VectorOfVariables

include("types.jl")

function Base.show(io::IO, opts::SolverOptions)
    longest_field_name = maximum([length(string(fname)) for fname in fieldnames(SolverOptions)])+2
    for name in fieldnames(SolverOptions)
        sname = string(name)
        pname = sname*repeat(" ", longest_field_name-length(sname))
        if getfield(opts,name) === nothing
            println(io, pname, ": NA")
        else
            println(io, pname, ": ", getfield(opts,name))
        end
    end
end

include("util.jl")
include("filter.jl")
include("printing.jl")
include("solver.jl")
include("model.jl")
include("BnBTree.jl")
include("MOI_wrapper/MOI_wrapper.jl")

end

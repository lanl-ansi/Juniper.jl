module Juniper

using MathProgBase
using JuMP
using JSON

include("types.jl")

function Base.show(io::IO, opts::SolverOptions) 
    longest_field_name = maximum([length(string(fname)) for fname in fieldnames(SolverOptions)])+2
    for name in fieldnames(SolverOptions)
        sname = string(name)
        pname = sname*repeat(" ", longest_field_name-length(sname))
        println(io, pname, ": ", getfield(opts,name))
    end
end

include("util.jl")
include("printing.jl")
include("solver.jl")
include("model.jl")
include("BnBTree.jl")


end

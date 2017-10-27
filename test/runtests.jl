using Logging

# suppress warnings during testing
Logging.configure(level=ERROR)

using Base.Test
using JuMP

using Ipopt
using PowerModels

using MINLPBnB
using MINLPBnB.BnBTree

include("load_mod.jl")

opt_rtol = 1e-6
opt_atol = 1e-6

sol_rtol = 1e-3
sol_atol = 1e-3

# model_from_mod("data/batch0812_nc.mod")

minlpbnb = MINLPBnB.MINLPBnBSolver(IpoptSolver(print_level=0);print_syms=[:NewIncumbent])

start = time()
include("basic.jl")
println("Time for all tests: ", time()-start)
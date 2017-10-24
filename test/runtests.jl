using Logging

reload("MINLPBnB")
reload("BnBTree")

# suppress warnings during testing
Logging.configure(level=ERROR)

using Base.Test
using JuMP

using Ipopt
using Gurobi

import MINLPBnB

opt_rtol = 1e-6
opt_atol = 1e-6

sol_rtol = 1e-3
sol_atol = 1e-3

gurobi = GurobiSolver()

minlpbnb = MINLPBnB.MINLPBnBSolver(IpoptSolver(print_level=0);print_syms=[:NewIncumbent])

start = time()
include("basic.jl")
println("Time for all tests: ", time()-start)
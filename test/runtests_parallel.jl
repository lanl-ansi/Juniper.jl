using Logging

# suppress warnings during testing
Logging.configure(level=ERROR)

using Base.Test
using JuMP

using Ipopt
using PowerModels

using MINLPBnB

opt_rtol = 1e-6
opt_atol = 1e-6

sol_rtol = 1e-3
sol_atol = 1e-3                           

start = time()

include("parallel.jl")
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

minlpbnb_strong_restart_2 = MINLPBnBSolver(IpoptSolver(print_level=0);
            branch_strategy=:StrongPseudoCost,
            strong_branching_nvars = 5,
            strong_branching_nsteps = 2,
            strong_restart = true
            )

minlpbnb_strong_restart = MINLPBnBSolver(IpoptSolver(print_level=0);
                branch_strategy=:StrongPseudoCost,
                strong_branching_nvars = 5,
                strong_restart = true
            )
minlpbnb_strong_no_restart = MINLPBnBSolver(IpoptSolver(print_level=0);
                branch_strategy=:StrongPseudoCost,
                strong_branching_nvars = 5,
                strong_restart = false
            )

minlpbnb_mosti = MINLPBnBSolver(IpoptSolver(print_level=0);
                branch_strategy=:MostInfeasible,
            )  

minlpbnb_pseudo = MINLPBnBSolver(IpoptSolver(print_level=0);
                branch_strategy=:PseudoCost,
            )                               

start = time()

# include("basic.jl")
include("user_limits.jl")
# include("pod.jl")
# include("power_models_acp.jl")
# include("power_models_socwr.jl")
println("Time for all tests: ", time()-start)
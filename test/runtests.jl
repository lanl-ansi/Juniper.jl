using Base,Logging

# suppress warnings during testing
Logging.configure(level=ERROR)

using Base.Test

if nworkers() > 1
    rmprocs(workers())
end
processors = 4 # Sys.CPU_CORES

if Base.JLOptions().code_coverage == 1
    addprocs(processors, exeflags = ["--code-coverage=user", "--inline=no", "--check-bounds=yes"])
else
    addprocs(processors, exeflags = "--check-bounds=yes")
end

println("Workers:", nworkers())

using JuMP

using Ipopt
using GLPKMathProgInterface
using PowerModels

using MINLPBnB

opt_rtol = 1e-6
opt_atol = 1e-6

sol_rtol = 1e-3
sol_atol = 1e-3

function DefaultTestSolver(;nl_solver=IpoptSolver(print_level=0), solver_args...)
    solver_args_dict = Dict{Symbol,Any}()
    solver_args_dict[:log_levels] = []
    for v in solver_args
        solver_args_dict[v[1]] = v[2]
    end
    return MINLPBnBSolver(nl_solver; solver_args_dict...)
end

minlpbnb_strong_restart_2 = DefaultTestSolver(
            branch_strategy=:StrongPseudoCost,
            strong_branching_perc = 25,
            strong_branching_nsteps = 2,
            strong_restart = true
            )

minlpbnb_strong_restart = DefaultTestSolver(
                branch_strategy=:StrongPseudoCost,
                strong_branching_perc = 25,
                strong_restart = true
            )
minlpbnb_strong_no_restart = DefaultTestSolver(
                branch_strategy=:StrongPseudoCost,
                strong_branching_perc = 25,
                strong_restart = false
            )

minlpbnb_mosti = DefaultTestSolver(
                branch_strategy=:MostInfeasible,
            )  

minlpbnb_pseudo = DefaultTestSolver(
                branch_strategy=:PseudoCost,
            )                               

start = time()

include("functions.jl")
include("basic.jl")
include("user_limits.jl")
include("parallel.jl")
include("fpump.jl")
include("pod.jl")
include("power_models_acp.jl")
include("power_models_socwr.jl")
println("Time for all tests: ", time()-start)
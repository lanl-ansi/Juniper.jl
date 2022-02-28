using Distributed

if nworkers() > 1
    rmprocs(workers())
end
processors = 4 # Sys.CPU_CORES

if Base.JLOptions().code_coverage == 1
    addprocs(
        processors,
        exeflags = [
            "--code-coverage=user",
            "--inline=no",
            "--check-bounds=yes",
        ],
    )
else
    addprocs(processors, exeflags = "--check-bounds=yes")
end

println("Workers:", nworkers())

using JuMP
using Test

import GLPK
import HiGHS
import Ipopt
import Juniper

const opt_rtol = 1e-6
const opt_atol = 1e-6
const sol_rtol = 1e-3
const sol_atol = 1e-3

function DefaultTestSolver(;
    nl_solver = optimizer_with_attributes(
        Ipopt.Optimizer,
        "print_level" => 0,
        "sb" => "yes",
    ),
    solver_args...,
)
    solver_args_result = Pair{String,Any}[]
    push!(solver_args_result, "log_levels" => Symbol[])
    push!(solver_args_result, "nl_solver" => nl_solver)
    for v in solver_args
        push!(solver_args_result, string(v[1]) => v[2])
    end
    return solver_args_result
end

const juniper_strong_restart_2 = DefaultTestSolver(
    branch_strategy = :StrongPseudoCost,
    strong_branching_perc = 25,
    strong_branching_nsteps = 2,
    strong_restart = true,
)

const juniper_reliable_restart = DefaultTestSolver(
    branch_strategy = :Reliability,
    reliability_branching_perc = 25,
    reliability_branching_threshold = 2,
    strong_restart = true,
)

const juniper_strong_restart = DefaultTestSolver(
    branch_strategy = :StrongPseudoCost,
    strong_branching_perc = 25,
    strong_restart = true,
)
const juniper_strong_no_restart = DefaultTestSolver(
    branch_strategy = :StrongPseudoCost,
    strong_branching_perc = 25,
    strong_restart = false,
)

const juniper_mosti = DefaultTestSolver(branch_strategy = :MostInfeasible)

const juniper_pseudo = DefaultTestSolver(branch_strategy = :PseudoCost)

start = time()

# TODO(odow): files not tested
#  * include("current_118.jl")
#  * include("power_models_acp.jl")
#  * include("power_models_socwr.jl")
#  * include("MINLPTests/run_minlptests")

@testset "Juniper" begin
    include("debug.jl")
    include("functions.jl")
    include("basic.jl")
    include("user_limits.jl")
    include("parallel.jl")
    include("fpump.jl")
    include("pod.jl")
    include("MOI_wrapper.jl")
end
println("Time for all tests: ", time() - start)

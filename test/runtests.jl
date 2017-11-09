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
using PowerModels

using MINLPBnB

opt_rtol = 1e-6
opt_atol = 1e-6

sol_rtol = 1e-3
sol_atol = 1e-3

function DefaultTestSolver(;nl_solver=IpoptSolver(print_level=0),
    log_levels                  = [],
    branch_strategy             = :StrongPseudoCost,
    # Strong branching
    strong_branching_nvars      = 5,
    strong_branching_nsteps     = 1,
    strong_restart              = true,
    # Obj cuts
    incumbent_constr            = true,
    # :UserLimit
    time_limit                  = Inf,  
    mip_gap                     = 0,
    best_obj_stop               = NaN,
    solution_limit              = 0,
    all_solutions               = false,
    list_of_solutions           = false,
    # Parallel
    processors                  = 1,
    # Traversing
    traverse_strategy           = :BFS
)
    return MINLPBnBSolver(nl_solver;
        log_levels = log_levels,
        branch_strategy = branch_strategy,
        strong_branching_nvars = strong_branching_nvars,
        strong_branching_nsteps = strong_branching_nsteps,
        strong_restart = strong_restart,
        incumbent_constr = incumbent_constr,
        time_limit = time_limit,
        mip_gap = mip_gap,
        best_obj_stop = best_obj_stop,
        solution_limit = solution_limit,
        all_solutions = all_solutions,
        list_of_solutions = list_of_solutions,
        processors = processors,
        traverse_strategy = traverse_strategy)
end

minlpbnb_strong_restart_2 = DefaultTestSolver(
            branch_strategy=:StrongPseudoCost,
            strong_branching_nvars = 5,
            strong_branching_nsteps = 2,
            strong_restart = true
            )

minlpbnb_strong_restart = DefaultTestSolver(
                branch_strategy=:StrongPseudoCost,
                strong_branching_nvars = 5,
                strong_restart = true
            )
minlpbnb_strong_no_restart = DefaultTestSolver(
                branch_strategy=:StrongPseudoCost,
                strong_branching_nvars = 5,
                strong_restart = false
            )

minlpbnb_mosti = DefaultTestSolver(
                branch_strategy=:MostInfeasible,
            )  

minlpbnb_pseudo = DefaultTestSolver(
                branch_strategy=:PseudoCost,
            )                               

start = time()

include("basic.jl")
include("user_limits.jl")
include("parallel.jl")
include("pod.jl")
include("power_models_acp.jl")
include("power_models_socwr.jl")
println("Time for all tests: ", time()-start)
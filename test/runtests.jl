using Base, Logging

using Test, Distributed




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



using LinearAlgebra
using Statistics
using Random


using JuMP

using Ipopt
using Cbc
using GLPK
# using PowerModels

using MathOptInterface 

const MOI = MathOptInterface 
const MOIU = MOI.Utilities

using Juniper
include("util.jl")

opt_rtol = 1e-6
opt_atol = 1e-6

sol_rtol = 1e-3
sol_atol = 1e-3



function solve(m::Model)
    JuMP.optimize!(m)
    bm = JuMP.backend(m)
    return MOI.get(bm, MOI.TerminationStatus()) 
end

function getsolvetime(m::Model)
    bm = JuMP.backend(m)
    return MOI.get(bm, MOI.SolveTime()) 
end

function internalmodel(m::Model)
    bm = JuMP.backend(m)
    return bm.optimizer.model.inner
end

function getobjgap(m::Model)
    bm = JuMP.backend(m)
    return MOI.get(bm, MOI.RelativeGap()) 
end

juniper_strong_restart_2 = DefaultTestSolver(
            branch_strategy=:StrongPseudoCost,
            strong_branching_perc = 25,
            strong_branching_nsteps = 2,
            strong_restart = true
            )

juniper_reliable_restart = DefaultTestSolver(
            branch_strategy=:Reliability,
            reliability_branching_perc = 25,
            reliability_branching_threshold = 2,
            strong_restart = true
            )

juniper_strong_restart = DefaultTestSolver(
                branch_strategy=:StrongPseudoCost,
                strong_branching_perc = 25,
                strong_restart = true
            )
juniper_strong_no_restart = DefaultTestSolver(
                branch_strategy=:StrongPseudoCost,
                strong_branching_perc = 25,
                strong_restart = false
            )

juniper_mosti = DefaultTestSolver(
                branch_strategy=:MostInfeasible,
            )

juniper_pseudo = DefaultTestSolver(
                branch_strategy=:PseudoCost,
            )

start = time()

@testset "Juniper" begin
    include("debug.jl")
    include("functions.jl")
    include("basic.jl")
    include("user_limits.jl")
    include("parallel.jl")
    include("fpump.jl")
    include("pod.jl")
    include("MOI_wrapper.jl")
    # include("power_models_acp.jl")
    # include("power_models_socwr.jl")
end
println("Time for all tests: ", time()-start)
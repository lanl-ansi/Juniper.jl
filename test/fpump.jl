include("POD_experiment/blend029.jl")
include("basic/gamsworld.jl")

@testset "fp tests" begin

@testset "FP: blend029" begin
    println("==================================")
    println("FP: blend029")
    println("==================================")

    m,objval = get_blend029()

    setsolver(m, DefaultTestSolver(
            branch_strategy=:MostInfeasible,
            strong_restart = true,
            feasibility_pump = true,
            time_limit = 5,
            mip_solver=GLPKSolverMIP()
    ))
    status = solve(m)

    @test MINLPBnB.getnsolutions(internalmodel(m)) >= 1
end

@testset "FP: cvxnonsep_nsig20r_problem" begin
    println("==================================")
    println("FP: cvxnonsep_nsig20r_problem")
    println("==================================")

    m,objval = get_blend029()

    setsolver(m, DefaultTestSolver(
            branch_strategy=:MostInfeasible,
            strong_restart = true,
            feasibility_pump = true,
            time_limit = 1,
            mip_solver=GLPKSolverMIP()
    ))
    status = solve(m)

    @test MINLPBnB.getnsolutions(internalmodel(m)) >= 1
end

end
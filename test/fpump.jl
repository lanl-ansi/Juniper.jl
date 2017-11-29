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

@testset "FP: no linear" begin
    println("==================================")
    println("FP: no linear")
    println("==================================")

    m = Model()
    @variable(m, x[1:10], Bin)
    @NLconstraint(m, x[1]^2+x[2]^2 == 0)
    @objective(m, Max, sum(x))

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
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
            feasibility_pump = true,
            time_limit = 5,
            mip_solver=GLPKSolverMIP()
    ))
    status = solve(m)

    @test Juniper.getnsolutions(internalmodel(m)) >= 1
end

@testset "FP: cvxnonsep_nsig20r_problem" begin
    println("==================================")
    println("FP: cvxnonsep_nsig20r_problem")
    println("==================================")

    m,objval = get_blend029()

    setsolver(m, DefaultTestSolver(
            branch_strategy=:MostInfeasible,
            feasibility_pump = true,
            time_limit = 1,
            mip_solver=GLPKSolverMIP()
    ))
    status = solve(m)

    @test Juniper.getnsolutions(internalmodel(m)) >= 1
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
            feasibility_pump = true,
            time_limit = 1,
            mip_solver=GLPKSolverMIP()
    ))
    status = solve(m)

    @test Juniper.getnsolutions(internalmodel(m)) >= 1
end

@testset "FP: Integer test" begin
    println("==================================")
    println("FP: Integer Test")
    println("==================================")
    m = Model()

    @variable(m, 1 <= x[1:4] <= 5, Int)

    @objective(m, Min, x[1])

    @constraint(m, x[1] >= 0.9)
    @constraint(m, x[1] <= 1.1)
    @NLconstraint(m, (x[1]-x[2])^2 >= 0.1)
    @NLconstraint(m, (x[2]-x[3])^2 >= 0.1)
    @NLconstraint(m, (x[1]-x[3])^2 >= 0.1)
    @NLconstraint(m, (x[1]-x[4])^2 >= 0.1)
    @NLconstraint(m, (x[2]-x[4])^2 >= 0.1)
    @NLconstraint(m, (x[3]-x[4])^2 >= 0.1)

    setsolver(m, DefaultTestSolver(
        branch_strategy=:MostInfeasible,
        feasibility_pump = true,
        time_limit = 1,
        mip_solver=GLPKSolverMIP()
    ))

    status = solve(m)
    @test Juniper.getnsolutions(internalmodel(m)) >= 1
end

@testset "FP: Integer test2" begin
    println("==================================")
    println("FP: Integer Test 2")
    println("==================================")
    m = Model()

    @variable(m, 0 <= x <= 10, Int)
    @variable(m, y >= 0)
    @variable(m, 0 <= u <= 10, Int)
    @variable(m, w == 1)

    @objective(m, Min, -3x - y)

    @constraint(m, 3x + 10 <= 20)
    @NLconstraint(m, y^2 <= u*w)
    @NLconstraint(m, x^2 >= u*w)

    setsolver(m, DefaultTestSolver(
        branch_strategy=:MostInfeasible,
        feasibility_pump = true,
        time_limit = 1,
        mip_solver=GLPKSolverMIP()
    ))

    status = solve(m)
    @test Juniper.getnsolutions(internalmodel(m)) >= 1
end

end
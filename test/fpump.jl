include("POD_experiment/tspn05.jl")
#include("POD_experiment/ndcc12persp.jl")
include("POD_experiment/FLay02H.jl")
include("basic/gamsworld.jl")

@testset "fp tests" begin

@testset "FP: no linear" begin
    println("==================================")
    println("FP: no linear")
    println("==================================")

    m = Model()
    @variable(m, x[1:10], Bin)
    @NLconstraint(m, x[1]^2+x[2]^2 == 0)
    @objective(m, Max, sum(x))

    set_optimizer(m, with_optimizer(
        Juniper.Optimizer,
        DefaultTestSolver(
            branch_strategy=:MostInfeasible,
            feasibility_pump = true,
            time_limit = 1,
            mip_solver=with_optimizer(Cbc.Optimizer, logLevel=0))
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

    set_optimizer(m, with_optimizer(
        Juniper.Optimizer,
        DefaultTestSolver(
            branch_strategy=:MostInfeasible,
            feasibility_pump = true,
            time_limit = 1,
            mip_solver=with_optimizer(Cbc.Optimizer, logLevel=0))
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

    set_optimizer(m, with_optimizer(
        Juniper.Optimizer,
        DefaultTestSolver(
            branch_strategy=:MostInfeasible,
            feasibility_pump = true,
            time_limit = 1,
            mip_solver=with_optimizer(Cbc.Optimizer, logLevel=0))
    ))

    status = solve(m)
    @test Juniper.getnsolutions(internalmodel(m)) >= 1
end

@testset "FP: infeasible cos" begin
    println("==================================")
    println("FP: Infeasible cos")
    println("==================================")
    m = Model()

    @variable(m, 1 <= x <= 5, Int)
    @variable(m, -2 <= y <= 2, Int)

    @objective(m, Min, -x-y)

    @NLconstraint(m, y==2*cos(2*x))

    set_optimizer(m, with_optimizer(
        Juniper.Optimizer,
        DefaultTestSolver(
            branch_strategy=:MostInfeasible,
            feasibility_pump = true,
            time_limit = 1,
            mip_solver=with_optimizer(Cbc.Optimizer, logLevel=0))
    ))

    status = solve(m)
    println("Status: ", status)

    @test status == MOI.LOCALLY_INFEASIBLE
    @test Juniper.getnsolutions(internalmodel(m)) == 0
end

@testset "FP: tspn05" begin
    println("==================================")
    println("FP: tspn05")
    println("==================================")

    m = get_tspn05()

    set_optimizer(m, with_optimizer(
        Juniper.Optimizer,
        DefaultTestSolver(
            branch_strategy=:StrongPseudoCost,
            feasibility_pump = true,
            mip_solver=with_optimizer(Cbc.Optimizer, logLevel=0))
    ))

    status = solve(m)

    @test status == MOI.LOCALLY_SOLVED
    @test isapprox(JuMP.objective_value(m),191.2541,atol=1e0)
end

@testset "FP: FLay02H" begin
    println("==================================")
    println("FP: FLay02H")
    println("==================================")

    # This probably needs a restart in real nlp
    m = get_FLay02H()

    set_optimizer(m, with_optimizer(
        Juniper.Optimizer,
        DefaultTestSolver(
            branch_strategy=:StrongPseudoCost,
            feasibility_pump = true,
            feasibility_pump_time_limit = 10,
            time_limit = 10,
            mip_solver=with_optimizer(Cbc.Optimizer, logLevel=0))
    ))

    status = solve(m)

    @test status == MOI.LOCALLY_SOLVED || status == MOI.TIME_LIMIT
    @test Juniper.getnsolutions(internalmodel(m)) >= 1
end

@testset "FP: FLay02H short feasibility_pump_time_limit" begin
    println("==================================")
    println("FP: FLay02H short feasibility_pump_time_limit")
    println("==================================")

    # This probably needs a restart in real nlp
    m = get_FLay02H()

    set_optimizer(m, with_optimizer(
        Juniper.Optimizer,
        DefaultTestSolver(
            branch_strategy=:StrongPseudoCost,
            feasibility_pump = true,
            feasibility_pump_time_limit = 1,
            time_limit = 2,
            mip_solver=with_optimizer(Cbc.Optimizer, logLevel=0))
    ))

    status = solve(m)

    @test status == MOI.LOCALLY_SOLVED || status == MOI.TIME_LIMIT
end

end
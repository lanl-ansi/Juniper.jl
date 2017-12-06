include("basic/gamsworld.jl")
include("POD_experiment/blend029.jl")

@testset "parallel tests" begin

@testset "Batch.mod reliable parallel" begin
    println("==================================")
    println("BATCH.MOD reliable")
    println("==================================")

    m = batch_problem()

    juniper = DefaultTestSolver(
        branch_strategy=:Reliability,
        strong_restart = false,
        processors = 4
    ) 

    setsolver(m, juniper)
    status = solve(m)
    @test status == :Optimal

    juniper_val = getobjectivevalue(m)
    juniper_bb = getobjbound(m)

    println("Solution by Juniper")
    println("obj: ", juniper_val)
    println("bound: ", juniper_bb)


    @test isapprox(juniper_val, 285506.5082, atol=opt_atol, rtol=opt_rtol)
end

@testset "Batch.mod reliable parallel > processors" begin
    println("==================================")
    println("BATCH.MOD reliable more processors than available")
    println("==================================")

    m = batch_problem()

    juniper = DefaultTestSolver(
        branch_strategy=:Reliability,
        strong_restart = false,
        processors = 10
    ) 

    setsolver(m, juniper)
    status = solve(m)
    @test status == :Optimal

    juniper_val = getobjectivevalue(m)
    juniper_bb = getobjbound(m)

    println("Solution by Juniper")
    println("obj: ", juniper_val)
    println("bound: ", juniper_bb)


    @test isapprox(juniper_val, 285506.5082, atol=opt_atol, rtol=opt_rtol)
end

@testset "Batch.mod no restart parallel" begin
    println("==================================")
    println("BATCH.MOD NO RESTART")
    println("==================================")

    m = batch_problem()

    juniper = DefaultTestSolver(
        branch_strategy=:StrongPseudoCost,
        strong_restart = false,
        processors = 4
    ) 

    setsolver(m, juniper)
    status = solve(m)
    @test status == :Optimal

    juniper_val = getobjectivevalue(m)
    juniper_bb = getobjbound(m)

    println("Solution by Juniper")
    println("obj: ", juniper_val)
    println("bound: ", juniper_bb)


    @test isapprox(juniper_val, 285506.5082, atol=opt_atol, rtol=opt_rtol)
end

@testset "Knapsack 100% limit" begin
    println("==================================")
    println("KNAPSACK 100%")
    println("==================================")

    m = Model(solver=DefaultTestSolver(;processors=2,traverse_strategy=:DBFS,mip_gap=100,
              branch_strategy=:MostInfeasible))
    
    v = [10,20,12,23,42]
    w = [12,45,12,22,21]
    @variable(m, x[1:5], Bin)

    @objective(m, Max, dot(v,x))

    @NLconstraint(m, sum(w[i]*x[i]^2 for i=1:5) <= 45)   

    status = solve(m)
    objval = getobjectivevalue(m)
    println("Obj: ", objval)
    best_bound_val = getobjbound(m)
    gap_val = getobjgap(m)
    println("bb: ", getobjbound(m))

    @test status == :UserLimit

    @test best_bound_val >= objval
    @test 0.01 <= gap_val <= 1 || Juniper.getnsolutions(internalmodel(m)) == 1
end

@testset "blend029" begin
    println("==================================")
    println("blend029")
    println("==================================")

    m,objval = get_blend029()

    setsolver(m, DefaultTestSolver(
            branch_strategy=:StrongPseudoCost,
            strong_branching_perc = 50,
            strong_branching_nsteps = 5,
            strong_restart = true,
            processors = 4
    ))
    status = solve(m)

    @test status == :Optimal

    juniper_val = getobjectivevalue(m)
    best_bound_val = getobjbound(m)
    gap_val = getobjgap(m)

    println("Solution by Juniper")
    println("obj: ", juniper_val)
    println("best_bound_val: ", best_bound_val)
    println("gap_val: ", gap_val)

    @test isapprox(juniper_val, objval, atol=1e0)
    @test isapprox(best_bound_val, objval, atol=1e0)
    @test isapprox(gap_val, 0, atol=1e-2)
end

@testset "bruteforce" begin
    println("==================================")
    println("Bruteforce")
    println("==================================")
    juniper_all_solutions = DefaultTestSolver(
        branch_strategy=:StrongPseudoCost,
        all_solutions = true,
        list_of_solutions = true,
        strong_restart = true,
        processors = 3
    )

    m = Model(solver=juniper_all_solutions)

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

    status = solve(m)
    println("Status: ", status)
    list_of_solutions = Juniper.getsolutions(internalmodel(m))
    @test length(unique(list_of_solutions)) == Juniper.getnsolutions(internalmodel(m))

    @test status == :Optimal
    @test Juniper.getnsolutions(internalmodel(m)) == 24
end

@testset "bruteforce 2 vs 1" begin
    println("==================================")
    println("Bruteforce 2 vs 1")
    println("==================================")
    juniper_all_solutions = DefaultTestSolver(
        branch_strategy=:PseudoCost,
        all_solutions = true,
        list_of_solutions = true,
        strong_restart = false,
        processors = 1
    )

    juniper_all_solutions_p2 = DefaultTestSolver(
        branch_strategy=:PseudoCost,
        all_solutions = true,
        list_of_solutions = true,
        strong_restart = false,
        processors = 2
    )

    m = Model(solver=juniper_all_solutions)

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

    status = solve(m)
    nbranches = Juniper.getnbranches(internalmodel(m))

    setsolver(m, juniper_all_solutions_p2)

    status = solve(m)
    @test Juniper.getnbranches(internalmodel(m)) == nbranches
end


@testset "bruteforce PseudoCost" begin
    println("==================================")
    println("Bruteforce PseudoCost")
    println("==================================")
    juniper_all_solutions = DefaultTestSolver(
        branch_strategy=:PseudoCost,
        all_solutions = true,
        list_of_solutions = true,
        strong_restart = true,
        processors = 3
    )

    m = Model(solver=juniper_all_solutions)

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

    status = solve(m)
    println("Status: ", status)
    list_of_solutions = Juniper.getsolutions(internalmodel(m))
    @test length(unique(list_of_solutions)) == Juniper.getnsolutions(internalmodel(m))

    @test status == :Optimal
    @test Juniper.getnsolutions(internalmodel(m)) == 24
end


@testset "time limit 5s" begin
    println("==================================")
    println("time imit 5s")
    println("==================================")
    juniper_all_solutions = DefaultTestSolver(
        branch_strategy=:PseudoCost,
        strong_restart = true,
        processors = 3,
        time_limit = 5
    )

    m,objval = get_blend029()
    setsolver(m, juniper_all_solutions)

    status = solve(m)
    println("Status: ", status)
    
    @test status == :UserLimit
    @test getsolvetime(m) <= 8 # some seconds more are allowed
end

@testset "infeasible cos" begin
    println("==================================")
    println("Infeasible cos")
    println("==================================")
    m = Model(solver=DefaultTestSolver(
        branch_strategy=:StrongPseudoCost,
        processors = 2
    ))

    @variable(m, 1 <= x <= 5, Int)
    @variable(m, -2 <= y <= 2, Int)

    @objective(m, Min, -x-y)

    @NLconstraint(m, y==2*cos(2*x))

    status = solve(m)
    println("Status: ", status)

    @test status == :Infeasible
end

@testset "infeasible relaxation" begin
    println("==================================")
    println("Infeasible relaxation")
    println("==================================")
    m = Model(solver=DefaultTestSolver(
        branch_strategy=:StrongPseudoCost,
        processors = 2
    ))

    @variable(m, 0 <= x[1:10] <= 2, Int)

    @objective(m, Min, sum(x))

    @constraint(m, sum(x[1:5]) <= 20)
    @NLconstraint(m, x[1]*x[2]*x[3] >= 10)

    status = solve(m)
    println("Status: ", status)

    @test status == :Infeasible
end

@testset "infeasible integer" begin
    println("==================================")
    println("Infeasible integer")
    println("==================================")
    m = Model(solver=DefaultTestSolver(
        branch_strategy=:StrongPseudoCost,
        processors = 2
    ))

    @variable(m, 0 <= x[1:10] <= 2, Int)

    @objective(m, Min, sum(x))

    @constraint(m, sum(x[1:5]) <= 20)
    @NLconstraint(m, x[1]*x[2]*x[3] >= 7)
    @NLconstraint(m, x[1]*x[2]*x[3] <= 7.5)

    status = solve(m)
    println("Status: ", status)

    @test status == :Infeasible
end

@testset "infeasible in strong" begin
    println("==================================")
    println("Infeasible in strong")
    println("==================================")
    m = Model(solver=DefaultTestSolver(
        branch_strategy=:StrongPseudoCost,
        processors = 2
    ))

    @variable(m, 0 <= x[1:5] <= 2, Int)

    @objective(m, Min, sum(x))

    @NLconstraint(m, x[3]^2 <= 2)
    @NLconstraint(m, x[3]^2 >= 1.2)

    status = solve(m)
    println("Status: ", status)

    @test status == :Infeasible
end

end
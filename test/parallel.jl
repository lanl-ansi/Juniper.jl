include("basic/gamsworld.jl")
include("POD_experiment/blend029.jl")

@testset "parallel tests" begin

@testset "Batch.mod no restart parallel" begin
    println("==================================")
    println("BATCH.MOD NO RESTART")
    println("==================================")

    m = batch_problem()

    minlpbnb = DefaultTestSolver(
        branch_strategy=:StrongPseudoCost,
        strong_restart = false,
        processors = 4
    ) 

    setsolver(m, minlpbnb)
    status = solve(m)
    @test status == :Optimal

    minlpbnb_val = getobjectivevalue(m)
    minlpbnb_bb = getobjbound(m)

    println("Solution by MINLPBnb")
    println("obj: ", minlpbnb_val)
    println("bound: ", minlpbnb_bb)


    @test isapprox(minlpbnb_val, 285506.5082, atol=opt_atol, rtol=opt_rtol)
end

@testset "blend029" begin
    println("==================================")
    println("blend029")
    println("==================================")

    m,objval = get_blend029()

    setsolver(m, DefaultTestSolver(
            branch_strategy=:StrongPseudoCost,
            strong_branching_nvars = 15,
            strong_branching_nsteps = 5,
            strong_restart = true,
            processors = 4
    ))
    status = solve(m)

    @test status == :Optimal

    minlpbnb_val = getobjectivevalue(m)
    best_bound_val = getobjbound(m)
    gap_val = getobjgap(m)

    println("Solution by MINLPBnb")
    println("obj: ", minlpbnb_val)
    println("best_bound_val: ", best_bound_val)
    println("gap_val: ", gap_val)

    @test isapprox(minlpbnb_val, objval, atol=1e0)
    @test isapprox(best_bound_val, objval, atol=1e0)
    @test isapprox(gap_val, 0, atol=1e-2)
end

@testset "bruteforce" begin
    println("==================================")
    println("Bruteforce")
    println("==================================")
    minlpbnb_all_solutions = DefaultTestSolver(
        branch_strategy=:StrongPseudoCost,
        all_solutions = true,
        list_of_solutions = true,
        strong_restart = true,
        processors = 3
    )

    m = Model(solver=minlpbnb_all_solutions)

    @variable(m, 1 <= x[1:4] <= 5, Int)


    @objective(m, Min, x[1])

    @constraint(m, x[1] == 1)
    @NLconstraint(m, (x[1]-x[2])^2 >= 0.1)
    @NLconstraint(m, (x[2]-x[3])^2 >= 0.1)
    @NLconstraint(m, (x[1]-x[3])^2 >= 0.1)
    @NLconstraint(m, (x[1]-x[4])^2 >= 0.1)
    @NLconstraint(m, (x[2]-x[4])^2 >= 0.1)
    @NLconstraint(m, (x[3]-x[4])^2 >= 0.1)

    status = solve(m)
    println("Status: ", status)
    list_of_solutions = MINLPBnB.getsolutions(internalmodel(m))
    @test length(unique(list_of_solutions)) == MINLPBnB.getnsolutions(internalmodel(m))

    @test status == :Optimal
    @test MINLPBnB.getnsolutions(internalmodel(m)) == 24
end

@testset "bruteforce PseudoCost" begin
    println("==================================")
    println("Bruteforce PseudoCost")
    println("==================================")
    minlpbnb_all_solutions = DefaultTestSolver(
        branch_strategy=:PseudoCost,
        all_solutions = true,
        list_of_solutions = true,
        strong_restart = true,
        processors = 3
    )

    m = Model(solver=minlpbnb_all_solutions)

    @variable(m, 1 <= x[1:4] <= 5, Int)


    @objective(m, Min, x[1])

    @constraint(m, x[1] == 1)
    @NLconstraint(m, (x[1]-x[2])^2 >= 0.1)
    @NLconstraint(m, (x[2]-x[3])^2 >= 0.1)
    @NLconstraint(m, (x[1]-x[3])^2 >= 0.1)
    @NLconstraint(m, (x[1]-x[4])^2 >= 0.1)
    @NLconstraint(m, (x[2]-x[4])^2 >= 0.1)
    @NLconstraint(m, (x[3]-x[4])^2 >= 0.1)

    status = solve(m)
    println("Status: ", status)
    list_of_solutions = MINLPBnB.getsolutions(internalmodel(m))
    @test length(unique(list_of_solutions)) == MINLPBnB.getnsolutions(internalmodel(m))

    @test status == :Optimal
    @test MINLPBnB.getnsolutions(internalmodel(m)) == 24
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
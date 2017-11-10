include("basic/gamsworld.jl")

@testset "basic tests" begin

@testset "bruteforce" begin
    println("==================================")
    println("Bruteforce")
    println("==================================")
    minlpbnb_all_solutions = DefaultTestSolver(
        branch_strategy=:StrongPseudoCost,
        all_solutions = true,
        list_of_solutions = true,
        strong_restart = true
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
        strong_restart = true
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
    m = Model(solver=minlpbnb_strong_restart)

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
    m = Model(solver=minlpbnb_strong_no_restart)

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
    m = Model(solver=minlpbnb_strong_no_restart)

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
    m = Model(solver=minlpbnb_strong_no_restart)

    @variable(m, 0 <= x[1:5] <= 2, Int)

    @objective(m, Min, sum(x))

    @NLconstraint(m, x[3]^2 <= 2)
    @NLconstraint(m, x[3]^2 >= 1.2)

    status = solve(m)
    println("Status: ", status)

    @test status == :Infeasible
end

@testset "One Integer small Strong" begin
    println("==================================")
    println("One Integer small Strong")
    println("==================================")
    m = Model(solver=minlpbnb_strong_no_restart)

    @variable(m, x >= 0, Int)
    @variable(m, y >= 0)
    @variable(m, 0 <= u <= 10, Int)
    @variable(m, w == 1)

    @objective(m, Min, -3x - y)

    @constraint(m, 3x + 10 <= 20)
    @NLconstraint(m, y^2 <= u*w)

    status = solve(m)
    println("Obj: ", getobjectivevalue(m))
    println("x: ", getvalue(x))
    println("y: ", getvalue(x))

    @test status == :Optimal
    @test isapprox(getobjectivevalue(m), -12.162277, atol=opt_atol)
    @test isapprox(getvalue(x), 3, atol=sol_atol)
    @test isapprox(getvalue(y), 3.162277, atol=sol_atol)
end

@testset "One Integer small MostInfeasible" begin
    println("==================================")
    println("One Integer small MostInfeasible")
    println("==================================")
    m = Model(solver=minlpbnb_mosti)

    @variable(m, x >= 0, Int)
    @variable(m, y >= 0)
    @variable(m, 0 <= u <= 10, Int)
    @variable(m, w == 1)

    @objective(m, Min, -3x - y)

    @constraint(m, 3x + 10 <= 20)
    @NLconstraint(m, y^2 <= u*w)

    status = solve(m)
    println("Obj: ", getobjectivevalue(m))
    println("x: ", getvalue(x))
    println("y: ", getvalue(x))

    @test status == :Optimal
    @test isapprox(getobjectivevalue(m), -12.162277, atol=opt_atol)
    @test isapprox(getvalue(x), 3, atol=sol_atol)
    @test isapprox(getvalue(y), 3.162277, atol=sol_atol)
end

@testset "One Integer small PseudoCost" begin
    println("==================================")
    println("One Integer small PseudoCost")
    println("==================================")
    m = Model(solver=minlpbnb_pseudo)

    @variable(m, x >= 0, Int)
    @variable(m, y >= 0)
    @variable(m, 0 <= u <= 10, Int)
    @variable(m, w == 1)

    @objective(m, Min, -3x - y)

    @constraint(m, 3x + 10 <= 20)
    @NLconstraint(m, y^2 <= u*w)

    status = solve(m)
    println("Obj: ", getobjectivevalue(m))
    println("x: ", getvalue(x))
    println("y: ", getvalue(x))

    @test status == :Optimal
    @test isapprox(getobjectivevalue(m), -12.162277, atol=opt_atol)
    @test isapprox(getvalue(x), 3, atol=sol_atol)
    @test isapprox(getvalue(y), 3.162277, atol=sol_atol)
end

@testset "Three Integers Small Strong" begin
    println("==================================")
    println("Three Integers Small Strong")
    println("==================================")
    m = Model(solver=minlpbnb_strong_no_restart)

    @variable(m, x >= 0, Int)
    @variable(m, y >= 0, Int)
    @variable(m, 0 <= u <= 10, Int)
    @variable(m, w == 1)

    @objective(m, Min, -3x - y)

    @constraint(m, 3x + 10 <= 20)
    @NLconstraint(m, y^2 <= u*w)

    status = solve(m)
    println("Obj: ", getobjectivevalue(m))

    @test status == :Optimal
    @test isapprox(getobjectivevalue(m), -12, atol=opt_atol)
    @test isapprox(getvalue(x), 3, atol=sol_atol)
    @test isapprox(getvalue(y), 3, atol=sol_atol)
end

@testset "Three Integers Small MostInfeasible" begin
    println("==================================")
    println("Three Integers Small MostInfeasible")
    println("==================================")

    m = Model(solver=minlpbnb_mosti)

    @variable(m, x >= 0, Int)
    @variable(m, y >= 0, Int)
    @variable(m, 0 <= u <= 10, Int)
    @variable(m, w == 1)

    @objective(m, Min, -3x - y)

    @constraint(m, 3x + 10 <= 20)
    @NLconstraint(m, y^2 <= u*w)

    status = solve(m)
    println("Obj: ", getobjectivevalue(m))

    @test status == :Optimal
    @test isapprox(getobjectivevalue(m), -12, atol=opt_atol)
    @test isapprox(getvalue(x), 3, atol=sol_atol)
    @test isapprox(getvalue(y), 3, atol=sol_atol)
end

@testset "Three Integers Small PseudoCost" begin
    println("==================================")
    println("Three Integers Small PseudoCost")
    println("==================================")

    m = Model(solver=minlpbnb_pseudo)

    @variable(m, x >= 0, Int)
    @variable(m, y >= 0, Int)
    @variable(m, 0 <= u <= 10, Int)
    @variable(m, w == 1)

    @objective(m, Min, -3x - y)

    @constraint(m, 3x + 10 <= 20)
    @NLconstraint(m, y^2 <= u*w)

    status = solve(m)
    println("Obj: ", getobjectivevalue(m))

    @test status == :Optimal
    @test isapprox(getobjectivevalue(m), -12, atol=opt_atol)
    @test isapprox(getvalue(x), 3, atol=sol_atol)
    @test isapprox(getvalue(y), 3, atol=sol_atol)
end

@testset "Knapsack Max" begin
    println("==================================")
    println("KNAPSACK")
    println("==================================")

    m = Model(solver=DefaultTestSolver(;traverse_strategy=:DBFS))

    v = [10,20,12,23,42]
    w = [12,45,12,22,21]
    @variable(m, x[1:5], Bin)

    @objective(m, Max, dot(v,x))
    
    @NLconstraint(m, sum(w[i]*x[i]^2 for i=1:5) <= 45)   

    status = solve(m)
    println("Obj: ", getobjectivevalue(m))

    @test status == :Optimal
    @test isapprox(getobjectivevalue(m), 65, atol=opt_atol)
    @test isapprox(getobjectivebound(m), 65, atol=opt_atol)
    @test isapprox(getvalue(x), [0,0,0,1,1], atol=sol_atol)
end

@testset "Integer at root" begin
    println("==================================")
    println("INTEGER AT ROOT")
    println("==================================")
    m = Model(solver=DefaultTestSolver())

    @variable(m, x[1:6] <= 1, Int)
    @constraint(m, x[1:6] .== 1)
    @objective(m, Max, sum(x))
    @NLconstraint(m, x[1]*x[2]*x[3]+x[4]*x[5]*x[6] <= 100)
    status = solve(m)
    @test status == :Optimal
    @test isapprox(getobjectivevalue(m), 6, atol=opt_atol)
    @test isapprox(getobjectivebound(m), 0, atol=opt_atol) # Ipopt return 0
end

@testset "Knapsack Max with epsilon" begin
    println("==================================")
    println("KNAPSACK with epsilon")
    println("==================================")

    m = Model(solver=DefaultTestSolver(;traverse_strategy=:DBFS,obj_epsilon=0.5))

    v = [10,20,12,23,42]
    w = [12,45,12,22,21]
    @variable(m, x[1:5], Bin)

    @objective(m, Max, dot(v,x))

    @NLconstraint(m, sum(w[i]*x[i]^2 for i=1:5) <= 45)   

    status = solve(m)
    println("Obj: ", getobjectivevalue(m))

    @test status == :Optimal
    @test isapprox(getobjectivevalue(m), 65, atol=opt_atol)
    @test isapprox(getvalue(x), [0,0,0,1,1], atol=sol_atol)
end

@testset "Knapsack Max with epsilon too strong" begin
    println("==================================")
    println("KNAPSACK with epsilon too strong")
    println("==================================")

    m = Model(solver=DefaultTestSolver(;traverse_strategy=:DBFS,obj_epsilon=0.1))

    v = [10,20,12,23,42]
    w = [12,45,12,22,21]
    @variable(m, x[1:5], Bin)

    @objective(m, Max, dot(v,x))

    @NLconstraint(m, sum(w[i]*x[i]^2 for i=1:5) <= 45)   

    status = solve(m)

    @test status == :Infeasible
end

@testset "Batch.mod Restart" begin
    println("==================================")
    println("BATCH.MOD RESTART")
    println("==================================")

    m = batch_problem()

    setsolver(m, minlpbnb_strong_restart)
    status = solve(m)
    @test status == :Optimal

    minlpbnb_val = getobjectivevalue(m)

    println("Solution by MINLPBnb")
    println("obj: ", minlpbnb_val)

    @test isapprox(minlpbnb_val, 285506.5082, atol=opt_atol, rtol=opt_rtol)
end

@testset "Batch.mod Restart 2 Levels" begin
    println("==================================")
    println("BATCH.MOD RESTART 2 LEVELS")
    println("==================================")

    m = batch_problem()

    setsolver(m, minlpbnb_strong_restart_2)
    status = solve(m)
    @test status == :Optimal

    minlpbnb_val = getobjectivevalue(m)

    println("Solution by MINLPBnb")
    println("obj: ", minlpbnb_val)

    @test isapprox(minlpbnb_val, 285506.5082, atol=opt_atol, rtol=opt_rtol)
end


@testset "cvxnonsep_nsig20r.mod restart" begin
    println("==================================")
    println("cvxnonsep_nsig20r.MOD RESTART")
    println("==================================")

    m = cvxnonsep_nsig20r_problem()

    setsolver(m, minlpbnb_strong_restart)
    status = solve(m)
    @test status == :Optimal

    minlpbnb_val = getobjectivevalue(m)

    println("Solution by MINLPBnb")
    println("obj: ", minlpbnb_val)

    @test isapprox(minlpbnb_val, 80.9493, atol=opt_atol, rtol=opt_rtol)
end

@testset "cvxnonsep_nsig20r.mod no restart" begin
    println("==================================")
    println("cvxnonsep_nsig20r.MOD NO RESTART")
    println("==================================")

    m = cvxnonsep_nsig20r_problem()

    setsolver(m, minlpbnb_strong_no_restart)
    status = solve(m)
    @test status == :Optimal

    minlpbnb_val = getobjectivevalue(m)

    println("Solution by MINLPBnb")
    println("obj: ", minlpbnb_val)

    @test isapprox(minlpbnb_val, 80.9493, atol=opt_atol, rtol=opt_rtol)
end

end
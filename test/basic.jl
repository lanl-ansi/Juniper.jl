include("basic/gamsworld.jl")

@testset "basic tests" begin

@testset "no objective" begin
    println("==================================")
    println("No objective")
    println("==================================")
    juniper = DefaultTestSolver(log_levels=[:Table])

    m = Model(solver=juniper)

    @variable(m, x, Int)

    @constraint(m, x >= 0)
    @constraint(m, x <= 5)
    @NLconstraint(m, x^2 >= 17)
   
    status = solve(m)
    @test status == :Optimal
    @test isapprox(getvalue(x), 5, atol=sol_atol)
    @test Juniper.getnsolutions(internalmodel(m)) == 1
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
        debug = true
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
    debugDict = internalmodel(m).debugDict
    @test getnstate(debugDict,:Integral) == 24
    @test different_hashes(debugDict) == true
    counter_test(debugDict,Juniper.getnbranches(internalmodel(m)))

    list_of_solutions = Juniper.getsolutions(internalmodel(m))
    @test length(unique(list_of_solutions)) == Juniper.getnsolutions(internalmodel(m))

    @test status == :Optimal
    @test Juniper.getnsolutions(internalmodel(m)) == 24
end

@testset "bruteforce full strong w/o restart" begin
    println("==================================")
    println("Bruteforce full strong w/o restart")
    println("==================================")
    juniper_all_solutions = DefaultTestSolver(
        branch_strategy=:StrongPseudoCost,
        all_solutions = true,
        list_of_solutions = true,
        strong_branching_perc = 100,
        strong_branching_nsteps = 100,
        strong_restart = false
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

    list_of_solutions = Juniper.getsolutions(internalmodel(m))
    @test length(unique(list_of_solutions)) == Juniper.getnsolutions(internalmodel(m))

    @test status == :Optimal
    @test Juniper.getnsolutions(internalmodel(m)) == 24
end

@testset "bruteforce approx time limit" begin
    println("==================================")
    println("Bruteforce  approx time limit")
    println("==================================")
    juniper_all_solutions = DefaultTestSolver(
        branch_strategy=:StrongPseudoCost,
        strong_branching_approx_time_limit=2,
        all_solutions = true,
        list_of_solutions = true,
        strong_restart = true
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

    list_of_solutions = Juniper.getsolutions(internalmodel(m))
    @test length(unique(list_of_solutions)) == Juniper.getnsolutions(internalmodel(m))

    @test status == :Optimal
    @test Juniper.getnsolutions(internalmodel(m)) == 24
end

@testset "bruteforce approx time limit reliable" begin
    println("==================================")
    println("Bruteforce  approx time reliable")
    println("==================================")
    juniper_all_solutions = DefaultTestSolver(
        branch_strategy=:Reliability,
        strong_branching_approx_time_limit=0.02,
        reliablility_branching_perc=100,
        all_solutions = true,
        list_of_solutions = true,
        strong_restart = true
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

    list_of_solutions = Juniper.getsolutions(internalmodel(m))
    @test length(unique(list_of_solutions)) == Juniper.getnsolutions(internalmodel(m))

    @test status == :Optimal
    @test Juniper.getnsolutions(internalmodel(m)) == 24
end

@testset "bruteforce PseudoCost" begin
    println("==================================")
    println("Bruteforce PseudoCost")
    println("==================================")
    juniper_all_solutions = DefaultTestSolver(
        branch_strategy=:PseudoCost,
        all_solutions = true,
        list_of_solutions = true,
        strong_restart = true
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

@testset "bruteforce Reliability" begin
    println("==================================")
    println("Bruteforce Reliability")
    println("==================================")
    juniper_all_solutions = DefaultTestSolver(
        branch_strategy=:Reliability,
        all_solutions = true,
        list_of_solutions = true,
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

@testset "no integer" begin
    println("==================================")
    println("no integer")
    println("==================================")
    m = Model(solver=juniper_strong_restart)

    println("Create variables/constr/obj")
    @variable(m, 1 <= x <= 5)
    @variable(m, -2 <= y <= 2)

    @objective(m, Min, -x-y)

    @NLconstraint(m, y==2*cos(2*x))
    println("before solve")
    status = solve(m)
    println("Status: ", status)

    @test status == :Optimal
end

@testset "infeasible cos" begin
    println("==================================")
    println("Infeasible cos")
    println("==================================")
    m = Model(solver=juniper_strong_restart)

    @variable(m, 1 <= x <= 5, Int)
    @variable(m, -2 <= y <= 2, Int)

    @objective(m, Min, -x-y)

    @NLconstraint(m, y==2*cos(2*x))

    status = solve(m)

    println("Status: ", status)

    @test status == :Infeasible
    @test isnan(getobjgap(m))
end


@testset "infeasible int reliable" begin
    println("==================================")
    println("Infeasible int reliable")
    println("==================================")
    m = Model(solver=juniper_reliable_restart)

    @variable(m, 1 <= x <= 5, Int)
    @variable(m, -2 <= y <= 2, Int)

    @objective(m, Min, -x-y)

    @NLconstraint(m, y >= sqrt(2))
    @NLconstraint(m, y <= sqrt(3))

    status = solve(m)
    println("Status: ", status)

    @test status == :Infeasible
    @test isnan(getobjgap(m))
end

@testset "infeasible sin with different bounds" begin
    println("==================================")
    println("Infeasible  sin with different bounds")
    println("==================================")
    m = Model()

    @variable(m, x <= 5, Int)
    @variable(m, y >= 2, Int)

    @objective(m, Min, -x-y)

    @NLconstraint(m, y==sin(x))

    setsolver(m, JuniperSolver(IpoptSolver(print_level=0),
        branch_strategy=:MostInfeasible,
        feasibility_pump = true,
        time_limit = 1,
        mip_solver=CbcSolver()
    ))
end

@testset "infeasible relaxation" begin
    println("==================================")
    println("Infeasible relaxation")
    println("==================================")
    m = Model(solver=DefaultTestSolver(;debug=true))

    @variable(m, 0 <= x[1:10] <= 2, Int)

    @objective(m, Min, sum(x))

    @constraint(m, sum(x[1:5]) <= 20)
    @NLconstraint(m, x[1]*x[2]*x[3] >= 10)

    status = solve(m)
    debug1 = m.internalModel.debugDict

    m = Model(solver=DefaultTestSolver(;debug=true))

    @variable(m, 0 <= x[1:10] <= 2, Int)

    @objective(m, Min, sum(x))

    @constraint(m, sum(x[1:5]) <= 20)
    @NLconstraint(m, x[1]*x[2]*x[3] >= 10)

    status = solve(m)

    debug2 = m.internalModel.debugDict
    opts = m.internalModel.options

    # should be deterministic
    @test debug1[:relaxation][:nrestarts] == debug2[:relaxation][:nrestarts] == opts.num_resolve_root_relaxation
    for i=1:debug1[:relaxation][:nrestarts]
        @test debug1[:relaxation][:restarts][i] == debug2[:relaxation][:restarts][i]
    end

    println("Status: ", status)

    @test status == :Infeasible
end

@testset "infeasible relaxation 2" begin
    println("==================================")
    println("Infeasible relaxation 2")
    println("==================================")
    m = Model(solver=juniper_strong_no_restart)

    @variable(m, x[1:3], Int)
    @variable(m, y)

    @objective(m, Max, sum(x))

    @NLconstraint(m, x[1]^2+x[2]^2+x[3]^2+y^2 <= 3)
    @NLconstraint(m, x[1]^2*x[2]^2*x[3]^2*y^2 >= 10)

    status = solve(m)
    println("Status: ", status)

    @test status == :Infeasible
end


@testset "infeasible integer" begin
    println("==================================")
    println("Infeasible integer")
    println("==================================")
    m = Model(solver=juniper_strong_no_restart)

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
    m = Model(solver=juniper_strong_no_restart)

    @variable(m, 0 <= x[1:5] <= 2, Int)

    @objective(m, Min, sum(x))

    @NLconstraint(m, x[3]^2 <= 2)
    @NLconstraint(m, x[3]^2 >= 1.2)

    status = solve(m)
    println("Status: ", status)

    @test status == :Infeasible
end

@testset "One Integer small Reliable" begin
    println("==================================")
    println("One Integer small Reliable")
    println("==================================")
    m = Model(solver=juniper_reliable_restart)

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

@testset "One Integer small Strong" begin
    println("==================================")
    println("One Integer small Strong")
    println("==================================")
    m = Model(solver=juniper_strong_no_restart)

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
    m = Model(solver=juniper_mosti)

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
    m = Model(solver=juniper_pseudo)

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
    m = Model(solver=juniper_strong_no_restart)

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

    m = Model(solver=juniper_mosti)

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

    m = Model(solver=juniper_pseudo)

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

    m = Model(solver=DefaultTestSolver(;traverse_strategy=:DBFS,
            incumbent_constr=true,mip_solver=CbcSolver(),
            strong_branching_approx_time_limit=1))

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

@testset "Knapsack Max Reliable" begin
    println("==================================")
    println("KNAPSACK Reliable no restart")
    println("==================================")

    m = Model(solver=DefaultTestSolver(;branch_strategy=:Reliability,
              strong_restart=false,strong_branching_approx_time_limit=1,gain_mu=0.5))
 
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

    setsolver(m, juniper_strong_restart)
    status = solve(m)
    @test status == :Optimal

    juniper_val = getobjectivevalue(m)

    println("Solution by Juniper")
    println("obj: ", juniper_val)

    @test isapprox(juniper_val, 285506.5082, atol=opt_atol, rtol=opt_rtol)
end

@testset "Batch.mod Restart 2 Levels" begin
    println("==================================")
    println("BATCH.MOD RESTART 2 LEVELS")
    println("==================================")

    m = batch_problem()

    setsolver(m, juniper_strong_restart_2)
    status = solve(m)
    @test status == :Optimal

    juniper_val = getobjectivevalue(m)

    println("Solution by Juniper")
    println("obj: ", juniper_val)

    @test isapprox(juniper_val, 285506.5082, atol=opt_atol, rtol=opt_rtol)
end


@testset "cvxnonsep_nsig20r.mod restart" begin
    println("==================================")
    println("cvxnonsep_nsig20r.MOD RESTART")
    println("==================================")

    m = cvxnonsep_nsig20r_problem()

    setsolver(m, juniper_strong_restart)
    status = solve(m)
    @test status == :Optimal

    juniper_val = getobjectivevalue(m)

    println("Solution by Juniper")
    println("obj: ", juniper_val)

    @test isapprox(juniper_val, 80.9493, atol=opt_atol, rtol=opt_rtol)
end

@testset "cvxnonsep_nsig20r.mod no restart" begin
    println("==================================")
    println("cvxnonsep_nsig20r.MOD NO RESTART")
    println("==================================")

    m = cvxnonsep_nsig20r_problem()

    setsolver(m, juniper_strong_no_restart)
    status = solve(m)
    @test status == :Optimal

    juniper_val = getobjectivevalue(m)

    println("Solution by Juniper")
    println("obj: ", juniper_val)

    @test isapprox(juniper_val, 80.9493, atol=opt_atol, rtol=opt_rtol)
end

@testset "Knapsack solution limit and table print test" begin
    println("==================================")
    println("Knapsack solution limit and table print test")
    println("==================================")

    juniper_one_solution = DefaultTestSolver(
        log_levels=[:Table],
        branch_strategy=:MostInfeasible,
        solution_limit=1
    )
    
    m = Model()
    v = [10,20,12,23,42]
    w = [12,45,12,22,21]
    @variable(m, x[1:5], Bin)

    @objective(m, Max, dot(v,x))

    @NLconstraint(m, sum(w[i]*x[i]^2 for i=1:5) <= 45)   

    setsolver(m, juniper_one_solution)
    status = solve(m)
    @test status == :UserLimit

    # maybe 2 found at the same time
    @test Juniper.getnsolutions(internalmodel(m)) <= 2
    @test Juniper.getnsolutions(internalmodel(m)) >= 1
end

@testset "bruteforce obj_epsilon" begin
    println("==================================")
    println("Bruteforce obj_epsilon")
    println("==================================")
    juniper_obj_eps = DefaultTestSolver(
        branch_strategy=:StrongPseudoCost,
        obj_epsilon=0.4
    )

    m = Model(solver=juniper_obj_eps)

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

    # maybe 2 found at the same time
    @test status == :Optimal
end

@testset "bruteforce best_obj_stop not reachable" begin
    println("==================================")
    println("Bruteforce best_obj_stop not reachable")
    println("==================================")
    juniper_best_obj_stop = DefaultTestSolver(
        branch_strategy=:StrongPseudoCost,
        best_obj_stop=0.8
    )

    m = Model(solver=juniper_best_obj_stop)

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

    # not possible to reach but should be solved anyway
    @test status == :Optimal
end

@testset "bruteforce best_obj_stop reachable" begin
    println("==================================")
    println("Bruteforce best_obj_stop reachable")
    println("==================================")
    juniper_one_solution = DefaultTestSolver(
        branch_strategy=:StrongPseudoCost,
        best_obj_stop=1
    )

    m = Model(solver=juniper_one_solution)

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

    # reachable and should break => UserLimit
    @test status == :UserLimit
    # maybe 2 found at the same time
    @test Juniper.getnsolutions(internalmodel(m)) <= 2
    @test Juniper.getnsolutions(internalmodel(m)) >= 1
end

end
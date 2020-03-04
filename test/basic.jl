include("basic/gamsworld.jl")

@testset "basic tests" begin

@testset "no objective and start value" begin
    println("==================================")
    println("No objective and start value")
    println("==================================")
    juniper = DefaultTestSolver(log_levels=[:Table])

    m = Model(optimizer_with_attributes(
        Juniper.Optimizer,
        juniper...)
    )

    @variable(m, x, Int, start=3)

    @constraint(m, x >= 0)
    @constraint(m, x <= 5)
    @NLconstraint(m, x^2 >= 17)
   
    status = solve(m)
    inner = internalmodel(m)
    @test JuMP.termination_status(m) == MOI.LOCALLY_SOLVED
    @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
    @test JuMP.dual_status(m) == MOI.FEASIBLE_POINT
    @test status == MOI.LOCALLY_SOLVED
    @test isapprox(JuMP.value(x), 5, atol=sol_atol)
    @test Juniper.getnsolutions(inner) == 1
    @test inner.primal_start[1] == 3
end


@testset "bruteforce" begin
    println("==================================")
    println("Bruteforce")
    println("==================================")

    m = Model()
    
    @variable(m, 1 <= x[1:4] <= 5, Int)

    special_minimizer_fct(x) = x
    grad(x) = 1.0
    grad2(x) = 0.0

    register_args = [:special_minimizer_fct, 1, special_minimizer_fct, grad, grad2]
    JuMP.register(m, register_args...)

    @NLobjective(m,Min,special_minimizer_fct(x[1]))

    juniper_all_solutions = DefaultTestSolver(
        branch_strategy=:StrongPseudoCost,
        list_of_solutions = true,
        strong_restart = true,
        debug = true,
        registered_functions=[Juniper.register(register_args...)]
    )

    JuMP.set_optimizer(m, optimizer_with_attributes(
        Juniper.Optimizer,
        juniper_all_solutions...
        )
    )

    @constraint(m, x[1] >= 0.9)
    @constraint(m, x[1] <= 1.1)
    @NLconstraint(m, (x[1]-x[2])^2 >= 0.1)
    @NLconstraint(m, (x[2]-x[3])^2 >= 0.1)
    @NLconstraint(m, (x[1]-x[3])^2 >= 0.1)
    @NLconstraint(m, (x[1]-x[4])^2 >= 0.1)
    @NLconstraint(m, (x[2]-x[4])^2 >= 0.1)
    @NLconstraint(m, (x[3]-x[4])^2 >= 0.1)

    JuMP.set_optimizer_attribute(m, "all_solutions", true)
    JuMP.optimize!(m)
    bm = JuMP.backend(m)
    status = MOI.get(bm, MOI.TerminationStatus()) 

    innermodel = bm.optimizer.model.inner
    debugDict = innermodel.debugDict
    @test getnstate(debugDict,:Integral) == 24
    @test different_hashes(debugDict) == true
    counter_test(debugDict,Juniper.getnbranches(innermodel))

    list_of_solutions = Juniper.getsolutions(innermodel)
    @test length(unique(list_of_solutions)) == Juniper.getnsolutions(innermodel)

    @test status == MOI.LOCALLY_SOLVED
    @test Juniper.getnsolutions(innermodel) == 24
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

    m = Model(optimizer_with_attributes(
        Juniper.Optimizer,
        juniper_all_solutions...)
    )

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

    @test status == MOI.LOCALLY_SOLVED
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

    m = Model(optimizer_with_attributes(
        Juniper.Optimizer,
        juniper_all_solutions...)
    )

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

    @test status == MOI.LOCALLY_SOLVED
    @test Juniper.getnsolutions(internalmodel(m)) == 24
end

@testset "bruteforce time limit reliable" begin
    println("==================================")
    println("Bruteforce  time reliable")
    println("==================================")
    juniper_all_solutions = DefaultTestSolver(
        branch_strategy=:Reliability,
        strong_branching_time_limit=0.02,
        reliability_branching_perc=100,
        all_solutions = true,
        list_of_solutions = true,
        strong_restart = true
    )

    optimizer = optimizer_with_attributes(Juniper.Optimizer, juniper_all_solutions...)
    
    m = Model(optimizer)
    moi_optimizer = JuMP.backend(m).optimizer.model
    
    # default => 1
    @test MOI.supports(moi_optimizer, MOI.NumberOfThreads()) == true
    MOI.set(moi_optimizer, MOI.NumberOfThreads(), nothing)

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
    @test MOI.get(moi_optimizer, MOI.NumberOfThreads()) == 1
    # all solutions are saved => nsolutions should equal length(solutions)
    @test MOI.get(moi_optimizer, MOI.ResultCount()) == Juniper.getnsolutions(internalmodel(m))

    @test status == MOI.LOCALLY_SOLVED
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

    m = Model(optimizer_with_attributes(
        Juniper.Optimizer,
        juniper_all_solutions...)
    )

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

    @test status == MOI.LOCALLY_SOLVED
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

    m = Model(optimizer_with_attributes(
        Juniper.Optimizer,
        juniper_all_solutions...)
    )

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

    @test status == MOI.LOCALLY_SOLVED
    @test Juniper.getnsolutions(internalmodel(m)) == 24
end

@testset "no integer" begin
    println("==================================")
    println("no integer")
    println("==================================")
    
    m = Model(optimizer_with_attributes(
        Juniper.Optimizer,
        juniper_strong_restart...)
    )

    println("Create variables/constr/obj")
    @variable(m, 1 <= x <= 5, start = 2.7)
    @variable(m, -2 <= y <= 2)

    @objective(m, Min, -x-y)

    @NLconstraint(m, y==2*cos(2*x))
    println("before solve")
    status = solve(m)
    println("Status: ", status)

    @test status == MOI.LOCALLY_SOLVED
end

@testset "infeasible cos" begin
    println("==================================")
    println("Infeasible cos")
    println("==================================")
    m = Model(optimizer_with_attributes(
        Juniper.Optimizer,
        juniper_strong_restart...)
    )

    @variable(m, 1 <= x <= 5, Int)
    @variable(m, -2 <= y <= 2, Int)

    @objective(m, Min, -x-y)

    @NLconstraint(m, y==2*cos(2*x))

    status = solve(m)

    println("Status: ", status)

    @test status == MOI.LOCALLY_INFEASIBLE
    @test JuMP.termination_status(m) == MOI.LOCALLY_INFEASIBLE
    @test JuMP.primal_status(m) == MOI.INFEASIBLE_POINT
    @test JuMP.dual_status(m) == MOI.INFEASIBLE_POINT
    @test isnan(getobjgap(m))
end


@testset "infeasible int reliable" begin
    println("==================================")
    println("Infeasible int reliable")
    println("==================================")
    m = Model(optimizer_with_attributes(
        Juniper.Optimizer,
        juniper_reliable_restart...)
    )

    @variable(m, 1 <= x <= 5, Int)
    @variable(m, -2 <= y <= 2, Int)

    @objective(m, Min, -x-y)

    @NLconstraint(m, y >= sqrt(2))
    @NLconstraint(m, y <= sqrt(3))

    status = solve(m)
    println("Status: ", status)

    @test status == MOI.LOCALLY_INFEASIBLE
    @test isnan(getobjgap(m))
end

@testset "infeasible sin with different bounds" begin
    println("==================================")
    println("Infeasible  sin with different bounds")
    println("==================================")
    m = Model(optimizer_with_attributes(
        Juniper.Optimizer,
        DefaultTestSolver(
            branch_strategy=:MostInfeasible,
            feasibility_pump = true,
            time_limit = 1,
            mip_solver=optimizer_with_attributes(Cbc.Optimizer)
        )...)
    )

    @variable(m, x, Int, start=3)
    @variable(m, y >= 2, Int)

    @constraint(m, 0 <= x <= 5)

    @objective(m, Min, -x-y)

    @NLconstraint(m, y==sin(x))

    status = solve(m)

    @test status == MOI.LOCALLY_INFEASIBLE
end

@testset "infeasible relaxation" begin
    println("==================================")
    println("Infeasible relaxation")
    println("==================================")
    
    m = Model(optimizer_with_attributes(
        Juniper.Optimizer,
        DefaultTestSolver(;debug=true)...)
    )

    @variable(m, 0 <= x[1:10] <= 2, Int)

    @objective(m, Min, sum(x))

    @constraint(m, sum(x[1:5]) <= 20)
    @NLconstraint(m, x[1]*x[2]*x[3] >= 10)

    status = solve(m)

    debug1 = internalmodel(m).debugDict

    m = Model(optimizer_with_attributes(
        Juniper.Optimizer,
        DefaultTestSolver(;debug=true)...)
    )

    @variable(m, 0 <= x[1:10] <= 2, Int)

    @objective(m, Min, sum(x))

    @constraint(m, sum(x[1:5]) <= 20)
    @NLconstraint(m, x[1]*x[2]*x[3] >= 10)

    status = solve(m)

    debug2 = internalmodel(m).debugDict
    opts = internalmodel(m).options

    # should be deterministic
    @test debug1[:relaxation][:nrestarts] == debug2[:relaxation][:nrestarts] == opts.num_resolve_root_relaxation
    for i=1:debug1[:relaxation][:nrestarts]
        @test debug1[:relaxation][:restarts][i] == debug2[:relaxation][:restarts][i]
    end

    println("Status: ", status)

    @test status == MOI.LOCALLY_INFEASIBLE
end

@testset "infeasible relaxation 2" begin
    println("==================================")
    println("Infeasible relaxation 2")
    println("==================================")
    
    m = Model(optimizer_with_attributes(
        Juniper.Optimizer,
        juniper_strong_no_restart...)
    )

    @variable(m, x[1:3], Int)
    @variable(m, y)

    @objective(m, Max, sum(x))

    @NLconstraint(m, x[1]^2+x[2]^2+x[3]^2+y^2 <= 3)
    @NLconstraint(m, x[1]^2*x[2]^2*x[3]^2*y^2 >= 10)

    status = solve(m)
    println("Status: ", status)

    @test status == MOI.LOCALLY_INFEASIBLE
end


@testset "infeasible integer" begin
    println("==================================")
    println("Infeasible integer")
    println("==================================")
    m = Model(optimizer_with_attributes(
        Juniper.Optimizer,
        juniper_strong_no_restart...)
    )

    @variable(m, 0 <= x[1:10] <= 2, Int)

    @objective(m, Min, sum(x))

    @constraint(m, sum(x[1:5]) <= 20)
    @NLconstraint(m, x[1]*x[2]*x[3] >= 7)
    @NLconstraint(m, x[1]*x[2]*x[3] <= 7.5)

    status = solve(m)
    println("Status: ", status)

    @test status == MOI.LOCALLY_INFEASIBLE
end

@testset "infeasible in strong" begin
    println("==================================")
    println("Infeasible in strong")
    println("==================================")
    m = Model(optimizer_with_attributes(
        Juniper.Optimizer,
        juniper_strong_no_restart...)
    )

    @variable(m, 0 <= x[1:5] <= 2, Int)

    @objective(m, Min, sum(x))

    @NLconstraint(m, x[3]^2 <= 2)
    @NLconstraint(m, x[3]^2 >= 1.2)

    status = solve(m)
    println("Status: ", status)

    @test status == MOI.LOCALLY_INFEASIBLE
end

@testset "One Integer small Reliable" begin
    println("==================================")
    println("One Integer small Reliable")
    println("==================================")
    m = Model()

    @variable(m, x >= 0, Int, start=2)
    @variable(m, y >= 0)
    @variable(m, 0 <= u <= 10, Int)
    @variable(m, w == 1)

    myf(x,y) = -3x-y
    function ∇f(g,x,y)
        g[1] = -3
        g[2] = -1
    end

    register_args = [:myf, 2, myf, ∇f]
    JuMP.register(m, register_args...)

    @NLobjective(m,Min,myf(x,y))

    juniper_reliable_restart_registered = DefaultTestSolver(
        branch_strategy=:Reliability,
        reliability_branching_perc = 25,
        reliability_branching_threshold = 2,
        strong_restart = true,
        registered_functions=[Juniper.register(register_args...)]
    )

    JuMP.set_optimizer(m, optimizer_with_attributes(
        Juniper.Optimizer,
        juniper_reliable_restart_registered...
        )
    )

    @constraint(m, 3x + 10 <= 20)
    @NLconstraint(m, y^2 <= u*w)

    status = solve(m)
    println("Obj: ", JuMP.objective_value(m))
    println("x: ", JuMP.value(x))
    println("y: ", JuMP.value(y))

    @test status == MOI.LOCALLY_SOLVED
    @test isapprox(JuMP.objective_value(m), -12.162277, atol=opt_atol)
    @test isapprox(JuMP.value(x), 3, atol=sol_atol)
    @test isapprox(JuMP.value(y), 3.162277, atol=sol_atol)
end

@testset "One Integer small Strong" begin
    println("==================================")
    println("One Integer small Strong")
    println("==================================")

    m = Model(optimizer_with_attributes(
        Juniper.Optimizer,
        juniper_strong_no_restart...)
    )

    @variable(m, x >= 0, Int)
    @variable(m, y >= 0)
    @variable(m, 0 <= u <= 10, Int)
    @variable(m, w == 1)

    @objective(m, Min, -3x - y)

    @constraint(m, 3x + 10 <= 20)
    @NLconstraint(m, y^2 <= u*w)

    status = solve(m)
    println("Obj: ", JuMP.objective_value(m))
    println("x: ", JuMP.value(x))
    println("y: ", JuMP.value(y))

    @test status == MOI.LOCALLY_SOLVED
    @test isapprox(JuMP.objective_value(m), -12.162277, atol=opt_atol)
    @test isapprox(JuMP.value(x), 3, atol=sol_atol)
    @test isapprox(JuMP.value(y), 3.162277, atol=sol_atol)
end

@testset "One Integer small MostInfeasible" begin
    println("==================================")
    println("One Integer small MostInfeasible")
    println("==================================")
    
    m = Model()    

    @variable(m, x >= 0, Int)
    @variable(m, y >= 0)
    @variable(m, 0 <= u <= 10, Int)
    @variable(m, w == 1)

    special_minimizer_fct(x,y) = -3x-y
    register_args = [:special_minimizer_fct, 2, special_minimizer_fct]
    JuMP.register(m, register_args...; autodiff=true)

    @NLobjective(m,Min,special_minimizer_fct(x, y))

    @constraint(m, 3x + 10 <= 20)
    @NLconstraint(m, y^2 <= u*w)

    JuMP.set_optimizer(m, optimizer_with_attributes(
        Juniper.Optimizer,
        DefaultTestSolver(
            branch_strategy=:MostInfeasible,
            registered_functions=[Juniper.register(register_args...; autodiff=true)]
        )...)
    )

    status = solve(m)
    println("Obj: ", JuMP.objective_value(m))
    println("x: ", JuMP.value(x))
    println("y: ", JuMP.value(y))

    @test status == MOI.LOCALLY_SOLVED
    @test isapprox(JuMP.objective_value(m), -12.162277, atol=opt_atol)
    @test isapprox(JuMP.value(x), 3, atol=sol_atol)
    @test isapprox(JuMP.value(y), 3.162277, atol=sol_atol)
end

@testset "One Integer small PseudoCost" begin
    println("==================================")
    println("One Integer small PseudoCost")
    println("==================================")

    m = Model(optimizer_with_attributes(
        Juniper.Optimizer,
        juniper_pseudo...)
    )

    @variable(m, x >= 0, Int)
    @variable(m, y >= 0)
    @variable(m, 0 <= u <= 10, Int)
    @variable(m, w == 1)

    @objective(m, Min, -3x - y)

    @constraint(m, 3x + 10 <= 20)
    @NLconstraint(m, y^2 <= u*w)

    status = solve(m)
    println("Obj: ", JuMP.objective_value(m))
    println("x: ", JuMP.value(x))
    println("y: ", JuMP.value(y))

    @test status == MOI.LOCALLY_SOLVED
    @test isapprox(JuMP.objective_value(m), -12.162277, atol=opt_atol)
    @test isapprox(JuMP.value(x), 3, atol=sol_atol)
    @test isapprox(JuMP.value(y), 3.162277, atol=sol_atol)
end

@testset "Three Integers Small Strong" begin
    println("==================================")
    println("Three Integers Small Strong")
    println("==================================")
    
    m = Model(optimizer_with_attributes(
        Juniper.Optimizer,
        juniper_strong_no_restart...)
    )

    @variable(m, x >= 0, Int)
    @variable(m, y >= 0, Int)
    @variable(m, 0 <= u <= 10, Int)
    @variable(m, w == 1)

    @objective(m, Min, -3x - y)

    @constraint(m, 3x + 10 <= 20)
    @NLconstraint(m, y^2 <= u*w)

    status = solve(m)
    println("Obj: ", JuMP.objective_value(m))

    @test status == MOI.LOCALLY_SOLVED
    @test isapprox(JuMP.objective_value(m), -12, atol=opt_atol)
    @test isapprox(JuMP.value(x), 3, atol=sol_atol)
    @test isapprox(JuMP.value(y), 3, atol=sol_atol)
end

@testset "Three Integers Small MostInfeasible" begin
    println("==================================")
    println("Three Integers Small MostInfeasible")
    println("==================================")

    m = Model(optimizer_with_attributes(
        Juniper.Optimizer,
        juniper_mosti...)
    )

    @variable(m, x >= 0, Int)
    @variable(m, y >= 0, Int)
    @variable(m, 0 <= u <= 10, Int)
    @variable(m, w == 1)

    @objective(m, Min, -3x - y)

    @constraint(m, 3x + 10 <= 20)
    @NLconstraint(m, y^2 <= u*w)

    status = solve(m)
    println("Obj: ", JuMP.objective_value(m))

    @test status == MOI.LOCALLY_SOLVED
    @test isapprox(JuMP.objective_value(m), -12, atol=opt_atol)
    @test isapprox(JuMP.value(x), 3, atol=sol_atol)
    @test isapprox(JuMP.value(y), 3, atol=sol_atol)
end

@testset "Three Integers Small PseudoCost" begin
    println("==================================")
    println("Three Integers Small PseudoCost")
    println("==================================")

    m = Model(optimizer_with_attributes(
        Juniper.Optimizer,
        juniper_pseudo...)
    )


    @variable(m, x >= 0, Int)
    @variable(m, y >= 0, Int)
    @variable(m, 0 <= u <= 10, Int)
    @variable(m, w == 1)

    @objective(m, Min, -3x - y)

    @constraint(m, 3x + 10 <= 20)
    @NLconstraint(m, y^2 <= u*w)

    status = solve(m)
    println("Obj: ", JuMP.objective_value(m))

    @test status == MOI.LOCALLY_SOLVED
    @test isapprox(JuMP.objective_value(m), -12, atol=opt_atol)
    @test isapprox(JuMP.value(x), 3, atol=sol_atol)
    @test isapprox(JuMP.value(y), 3, atol=sol_atol)
end

@testset "Knapsack Max" begin
    println("==================================")
    println("KNAPSACK")
    println("==================================")

    m = Model(optimizer_with_attributes(
        Juniper.Optimizer,
        DefaultTestSolver(;traverse_strategy=:DBFS,
            incumbent_constr=true,mip_solver=optimizer_with_attributes(Cbc.Optimizer),
            strong_branching_time_limit=1)...)
    )

    v = [10,20,12,23,42]
    w = [12,45,12,22,21]
    @variable(m, x[1:5], Bin)

    @objective(m, Max, dot(v,x))
    
    @NLconstraint(m, sum(w[i]*x[i]^2 for i=1:5) <= 45)   

    status = solve(m)
    println("Obj: ", JuMP.objective_value(m))

    @test status == MOI.LOCALLY_SOLVED
    @test isapprox(JuMP.objective_value(m), 65, atol=opt_atol)
    @test isapprox(JuMP.objective_bound(m), 65, atol=opt_atol)
    @test isapprox(JuMP.value.(x), [0,0,0,1,1], atol=sol_atol)
end


@testset "Knapsack Max Reliable" begin
    println("==================================")
    println("KNAPSACK Reliable no restart")
    println("==================================")
 
    m = Model(optimizer_with_attributes(
        Juniper.Optimizer,
        DefaultTestSolver(;branch_strategy=:Reliability,
              strong_restart=false,strong_branching_time_limit=1,gain_mu=0.5)...)
    )


    v = [10,20,12,23,42]
    w = [12,45,12,22,21]
    @variable(m, x[1:5], Bin)

    @objective(m, Max, dot(v,x))

    @NLconstraint(m, sum(w[i]*x[i]^2 for i=1:5) <= 45)   

    status = solve(m)
    println("Obj: ", JuMP.objective_value(m))

    @test status == MOI.LOCALLY_SOLVED
    @test isapprox(JuMP.objective_value(m), 65, atol=opt_atol)
    @test isapprox(JuMP.objective_bound(m), 65, atol=opt_atol)
    @test isapprox(JuMP.value.(x), [0,0,0,1,1], atol=sol_atol)
end


@testset "Integer at root" begin
    println("==================================")
    println("INTEGER AT ROOT")
    println("==================================")
    m = Model(optimizer_with_attributes(
        Juniper.Optimizer,
        DefaultTestSolver()...)
    )

    @variable(m, x[1:6] <= 1, Int)
    @constraint(m, x[1:6] .== 1)
    @objective(m, Max, sum(x))
    @NLconstraint(m, x[1]*x[2]*x[3]+x[4]*x[5]*x[6] <= 100)
    status = solve(m)
    @test status == MOI.LOCALLY_SOLVED
    @test isapprox(JuMP.objective_value(m), 6, atol=opt_atol)
    # objective bound doesn't exists anymore in Ipopt
    #@test isapprox(JuMP.objective_bound(m), 0, atol=opt_atol) # Ipopt return 0
end

@testset "Knapsack Max with epsilon" begin
    println("==================================")
    println("KNAPSACK with epsilon")
    println("==================================")

    m = Model(optimizer_with_attributes(
        Juniper.Optimizer,
        DefaultTestSolver(;traverse_strategy=:DBFS,obj_epsilon=0.5)...)
    )

    v = [10,20,12,23,42]
    w = [12,45,12,22,21]
    @variable(m, x[1:5], Bin)

    @objective(m, Max, dot(v,x))

    @NLconstraint(m, sum(w[i]*x[i]^2 for i=1:5) <= 45)   

    status = solve(m)
    println("Obj: ", JuMP.objective_value(m))

    @test status == MOI.LOCALLY_SOLVED
    @test isapprox(JuMP.objective_value(m), 65, atol=opt_atol)
    @test isapprox(JuMP.value.(x), [0,0,0,1,1], atol=sol_atol)
end

@testset "Knapsack Max with epsilon too strong" begin
    println("==================================")
    println("KNAPSACK with epsilon too strong")
    println("==================================")

    m = Model(optimizer_with_attributes(
        Juniper.Optimizer,
        DefaultTestSolver(;traverse_strategy=:DBFS,obj_epsilon=0.1)...)
    )

    v = [10,20,12,23,42]
    w = [12,45,12,22,21]
    # set all to 1 which is infeasible otherwise incumbent solution would be found
    @variable(m, x[1:5], Bin, start=1)

    @objective(m, Max, dot(v,x))

    @NLconstraint(m, sum(w[i]*x[i]^2 for i=1:5) <= 45)   

    status = solve(m)

    @test status == MOI.LOCALLY_INFEASIBLE
end

@testset "Batch.mod Restart" begin
    println("==================================")
    println("BATCH.MOD RESTART")
    println("==================================")

    m = batch_problem()

    JuMP.set_optimizer(m, optimizer_with_attributes(
        Juniper.Optimizer,
        juniper_strong_restart...)
    )
    
    status = solve(m)
    @test status == MOI.LOCALLY_SOLVED

    juniper_val = JuMP.objective_value(m)

    println("Solution by Juniper")
    println("obj: ", juniper_val)

    @test isapprox(juniper_val, 285506.5082, atol=opt_atol, rtol=opt_rtol)
end

@testset "Batch.mod Restart 2 Levels" begin
    println("==================================")
    println("BATCH.MOD RESTART 2 LEVELS")
    println("==================================")

    m = batch_problem()

    JuMP.set_optimizer(m, optimizer_with_attributes(
        Juniper.Optimizer,
        juniper_strong_restart_2...)
    )

    status = solve(m)
    @test status == MOI.LOCALLY_SOLVED

    juniper_val = JuMP.objective_value(m)

    println("Solution by Juniper")
    println("obj: ", juniper_val)

    @test isapprox(juniper_val, 285506.5082, atol=opt_atol, rtol=opt_rtol)
end


@testset "cvxnonsep_nsig20r.mod restart" begin
    println("==================================")
    println("cvxnonsep_nsig20r.MOD RESTART")
    println("==================================")

    m = cvxnonsep_nsig20r_problem()

    JuMP.set_optimizer(m, optimizer_with_attributes(
        Juniper.Optimizer,
        juniper_strong_restart...)
    )

    status = solve(m)
    @test status == MOI.LOCALLY_SOLVED

    juniper_val = JuMP.objective_value(m)

    println("Solution by Juniper")
    println("obj: ", juniper_val)

    @test isapprox(juniper_val, 80.9493, atol=opt_atol, rtol=opt_rtol)
end

@testset "cvxnonsep_nsig20r.mod no restart" begin
    println("==================================")
    println("cvxnonsep_nsig20r.MOD NO RESTART")
    println("==================================")

    m = cvxnonsep_nsig20r_problem()

    JuMP.set_optimizer(m, optimizer_with_attributes(
        Juniper.Optimizer,
        juniper_strong_no_restart...)
    )

    status = solve(m)
    @test status == MOI.LOCALLY_SOLVED

    juniper_val = JuMP.objective_value(m)

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

    JuMP.set_optimizer(m, optimizer_with_attributes(
        Juniper.Optimizer,
        juniper_one_solution...)
    )

    status = solve(m)
    @test status == MOI.SOLUTION_LIMIT

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


    m = Model(optimizer_with_attributes(
        Juniper.Optimizer,
        juniper_obj_eps...)
    )

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
    @test status == MOI.LOCALLY_SOLVED
end

@testset "bruteforce best_obj_stop not reachable" begin
    println("==================================")
    println("Bruteforce best_obj_stop not reachable")
    println("==================================")
    juniper_best_obj_stop = DefaultTestSolver(
        branch_strategy=:StrongPseudoCost,
        best_obj_stop=0.8
    )
    m = Model(optimizer_with_attributes(
        Juniper.Optimizer,
        juniper_best_obj_stop...)
    )

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
    @test status == MOI.LOCALLY_SOLVED
end

@testset "bruteforce best_obj_stop reachable" begin
    println("==================================")
    println("Bruteforce best_obj_stop reachable")
    println("==================================")
    juniper_one_solution = DefaultTestSolver(
        branch_strategy=:StrongPseudoCost,
        best_obj_stop=1
    )

    m = Model(optimizer_with_attributes(
        Juniper.Optimizer,
        juniper_one_solution...)
    )

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

    # reachable and should break
    @test status == MOI.OBJECTIVE_LIMIT
    # maybe 2 found at the same time
    @test Juniper.getnsolutions(internalmodel(m)) <= 2
    @test Juniper.getnsolutions(internalmodel(m)) >= 1
end

# this test has a lot "Only almost locally solved" warnings (in mumps at least)
@testset "Sum 1/x = 2" begin
    println("==================================")
    println("Sum 1/x = 2")
    println("==================================")
    m = Model(optimizer_with_attributes(
        Juniper.Optimizer,
        juniper_pseudo...)
    )

    @variable(m, 1 <= x[1:11], Int)
    @constraint(m, [i=2:11], x[i-1] <= x[i] - 1)
    @NLconstraint(m, sum(1 / x[i] for i in 1:11) == 2)

    status = solve(m)

    # in mumps/on travis this returns ALMOST_LOCALLY_SOLVED but with ma27/locally it is LOCALLY_SOLVED
    @test Juniper.state_is_optimal(status; allow_almost=true)
    @test isapprox(2, sum(1/v for v in JuMP.value.(x)), atol=opt_atol)
end

# this test has a lot "Only almost locally solved" warnings
@testset "Sum 1/x = 2 don't allow almost" begin
    println("==================================")
    println("Sum 1/x = 2 don't allow almost")
    println("==================================")
    m = Model(optimizer_with_attributes(
        Juniper.Optimizer,
        DefaultTestSolver(
                branch_strategy=:PseudoCost,
                allow_almost_solved = false
            )...)
    )

    @variable(m, 1 <= x[1:11], Int)
    @constraint(m, [i=2:11], x[i-1] <= x[i] - 1)
    @NLconstraint(m, sum(1 / x[i] for i in 1:11) == 2)

    status = solve(m)

    # even without almost this should still be solveable
    @test status == MOI.LOCALLY_SOLVED
    @test isapprox(2, sum(1/v for v in JuMP.value.(x)), atol=opt_atol)
end

#this test has an expression where a variable will be dereferenced twice
@testset "Nested variable reference" begin
    println("==================================")
    println("Nested variable reference")
    println("==================================")

    m = Model(optimizer_with_attributes(
        Juniper.Optimizer,
        DefaultTestSolver()...)
    )

    x = @variable(m, x[i=1:5,j=1:5], Bin)
    xrowsum = @NLexpression(m, xrowsum[i=1:5], sum(x[i,j] for j in 1:5))
    @NLobjective(m, Max, sum(sum(x[i,j] + xrowsum[i] for i=1:5) for j = 1:5))

    optimize!(m)

    @test JuMP.termination_status(m) == MOI.LOCALLY_SOLVED
end

@testset "Expr dereferencing for @NLexpression" begin
    # See issue 184
    m = Model()
    set_optimizer(m, optimizer_with_attributes(
        Juniper.Optimizer,
        DefaultTestSolver(
            mip_solver=optimizer_with_attributes(Cbc.Optimizer, "logLevel" => 0))...
    ))

    @variable(m, 0 <= a_var <= 1) 
    @variable(m, bin_var, Int) 
    b_expr = @NLexpression(m, bin_var/1)
    @NLconstraint(m, 1.1 >= b_expr) 
    an_expr = @NLexpression(m, a_var / 1) 

    @NLobjective(m, Max, bin_var + an_expr)

    optimize!(m)
    @test JuMP.objective_value(m) ≈ 2.0
    @test JuMP.value(bin_var) ≈ 1
    @test JuMP.value(a_var) ≈ 1
end

end
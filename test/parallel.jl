include("basic/gamsworld.jl")

@testset "parallel tests" begin
    @testset "Batch.mod reliable parallel" begin
        println("==================================")
        println("BATCH.MOD reliable")
        println("==================================")
        m = batch_problem()
        juniper = DefaultTestSolver(
            branch_strategy = :Reliability,
            strong_restart = false,
            processors = 4,
            mip_solver = optimizer_with_attributes(
                HiGHS.Optimizer,
                "output_flag" => false,
            ),
            incumbent_constr = true,
        )
        set_optimizer(
            m,
            optimizer_with_attributes(Juniper.Optimizer, juniper...),
        )
        optimize!(m)
        @test termination_status(m) == MOI.LOCALLY_SOLVED
        juniper_val = JuMP.objective_value(m)
        juniper_bb = JuMP.objective_bound(m)
        println("Solution by Juniper")
        println("obj: ", juniper_val)
        println("bound: ", juniper_bb)
        @test isapprox(
            juniper_val,
            285506.5082,
            atol = opt_atol,
            rtol = opt_rtol,
        )
    end

    @testset "Knapsack solution limit and table print test" begin
        println("==================================")
        println("Knapsack solution limit and table print test")
        println("==================================")
        juniper_one_solution = DefaultTestSolver(
            log_levels = [:Table],
            branch_strategy = :MostInfeasible,
            solution_limit = 1,
            mip_solver = optimizer_with_attributes(
                HiGHS.Optimizer,
                "output_flag" => false,
            ),
            processors = 3,
        )
        m = Model()
        v = [10, 20, 12, 23, 42]
        w = [12, 45, 12, 22, 21]
        @variable(m, x[1:5], Bin)
        @objective(m, Max, v' * x)
        @NLconstraint(m, sum(w[i] * x[i]^2 for i in 1:5) <= 45)
        set_optimizer(
            m,
            optimizer_with_attributes(
                Juniper.Optimizer,
                juniper_one_solution...,
            ),
        )
        optimize!(m)
        @test termination_status(m) == MOI.SOLUTION_LIMIT
        @test 1 <= result_count(m) <= 3
    end

    @testset "Knapsack Max Reliable incumbent_constr" begin
        println("==================================")
        println("KNAPSACK Reliable incumbent_constr")
        println("==================================")
        m = Model(
            optimizer_with_attributes(
                Juniper.Optimizer,
                DefaultTestSolver(;
                    branch_strategy = :MostInfeasible,
                    incumbent_constr = true,
                    processors = 2,
                )...,
            ),
        )
        v = [10, 20, 12, 23, 42]
        w = [12, 45, 12, 22, 21]
        @variable(m, x[1:5], Bin)
        @objective(m, Max, v' * x)
        @NLconstraint(m, sum(w[i] * x[i]^2 for i in 1:5) <= 45)
        optimize!(m)
        println("Obj: ", JuMP.objective_value(m))
        @test termination_status(m) == MOI.LOCALLY_SOLVED
        @test isapprox(JuMP.objective_value(m), 65, atol = opt_atol)
        @test isapprox(JuMP.value.(x), [0, 0, 0, 1, 1], atol = sol_atol)
    end

    @testset "Batch.mod reliable parallel > processors" begin
        println("==================================")
        println("BATCH.MOD reliable more processors than available")
        println("==================================")
        m = batch_problem()
        juniper = DefaultTestSolver(
            branch_strategy = :Reliability,
            strong_restart = false,
            processors = 10,
        )
        set_optimizer(
            m,
            optimizer_with_attributes(Juniper.Optimizer, juniper...),
        )
        optimize!(m)
        @test termination_status(m) == MOI.LOCALLY_SOLVED
        juniper_val = JuMP.objective_value(m)
        juniper_bb = JuMP.objective_bound(m)
        println("Solution by Juniper")
        println("obj: ", juniper_val)
        println("bound: ", juniper_bb)
        # must have changed to 4 processors
        @test unsafe_backend(m).inner.options.processors == 4
        @test isapprox(
            juniper_val,
            285506.5082,
            atol = opt_atol,
            rtol = opt_rtol,
        )
    end

    @testset "Batch.mod no restart parallel" begin
        println("==================================")
        println("BATCH.MOD NO RESTART")
        println("==================================")
        m = batch_problem()
        juniper = DefaultTestSolver(
            branch_strategy = :StrongPseudoCost,
            strong_restart = false,
            processors = 4,
        )
        set_optimizer(
            m,
            optimizer_with_attributes(Juniper.Optimizer, juniper...),
        )
        optimize!(m)
        @test termination_status(m) == MOI.LOCALLY_SOLVED
        juniper_val = JuMP.objective_value(m)
        juniper_bb = JuMP.objective_bound(m)
        println("Solution by Juniper")
        println("obj: ", juniper_val)
        println("bound: ", juniper_bb)
        @test isapprox(
            juniper_val,
            285506.5082,
            atol = opt_atol,
            rtol = opt_rtol,
        )
    end

    @testset "Knapsack 100% limit" begin
        println("==================================")
        println("KNAPSACK 100%")
        println("==================================")
        m = Model(
            optimizer_with_attributes(
                Juniper.Optimizer,
                DefaultTestSolver(;
                    processors = 2,
                    traverse_strategy = :DBFS,
                    mip_gap = 100,
                    branch_strategy = :MostInfeasible,
                )...,
            ),
        )
        v = [10, 20, 12, 23, 42]
        w = [12, 45, 12, 22, 21]
        # don't allow start to be an incumbent
        @variable(m, x[1:5], Bin, start = 0.5)
        @objective(m, Max, v' * x)
        @NLconstraint(m, sum(w[i] * x[i]^2 for i in 1:5) <= 45)
        optimize!(m)
        objval = JuMP.objective_value(m)
        println("Obj: ", objval)
        best_bound_val = JuMP.objective_bound(m)
        gap_val = relative_gap(m)
        println("bb: ", JuMP.objective_bound(m))
        @test termination_status(m) == MOI.OBJECTIVE_LIMIT
        @test best_bound_val >= objval
        @test 0.01 <= gap_val <= 1 || result_count(m) == 1
    end

    @testset "bruteforce" begin
        println("==================================")
        println("Bruteforce")
        println("==================================")
        juniper_all_solutions = DefaultTestSolver(
            branch_strategy = :StrongPseudoCost,
            all_solutions = true,
            list_of_solutions = true,
            strong_restart = true,
            processors = 3,
            debug = true,
        )
        m = Model(
            optimizer_with_attributes(
                Juniper.Optimizer,
                juniper_all_solutions...,
            ),
        )
        @variable(m, 1 <= x[1:4] <= 5, Int)
        @objective(m, Min, x[1])
        @constraint(m, x[1] >= 0.9)
        @constraint(m, x[1] <= 1.1)
        @NLconstraint(m, (x[1] - x[2])^2 >= 0.1)
        @NLconstraint(m, (x[2] - x[3])^2 >= 0.1)
        @NLconstraint(m, (x[1] - x[3])^2 >= 0.1)
        @NLconstraint(m, (x[1] - x[4])^2 >= 0.1)
        @NLconstraint(m, (x[2] - x[4])^2 >= 0.1)
        @NLconstraint(m, (x[3] - x[4])^2 >= 0.1)
        optimize!(m)
        debugDict = unsafe_backend(m).inner.debugDict
        list_of_solutions = Juniper.getsolutions(unsafe_backend(m).inner)
        @test length(unique(list_of_solutions)) == result_count(m)
        @test termination_status(m) == MOI.LOCALLY_SOLVED
        @test result_count(m) == 24
        @test getnstate(debugDict, :Integral) == 24
        @test different_hashes(debugDict) == true
        counter_test(debugDict, Juniper.getnbranches(unsafe_backend(m).inner))
    end

    @testset "bruteforce fake parallel vs sequential" begin
        println("==================================")
        println("Bruteforce fake parallel vs sequential")
        println("==================================")
        juniper_all_solutions = DefaultTestSolver(
            branch_strategy = :PseudoCost,
            all_solutions = true,
            list_of_solutions = true,
            strong_restart = false,
            processors = 1,
        )
        juniper_all_solutions_p = DefaultTestSolver(
            branch_strategy = :PseudoCost,
            all_solutions = true,
            list_of_solutions = true,
            strong_restart = false,
            processors = 1,
            force_parallel = true, # just for testing this goes into the parallel branch (using driver + 1)
        )
        m = Model(
            optimizer_with_attributes(
                Juniper.Optimizer,
                juniper_all_solutions...,
            ),
        )
        @variable(m, 1 <= x[1:4] <= 5, Int)
        @objective(m, Min, x[1])
        @constraint(m, x[1] >= 0.9)
        @constraint(m, x[1] <= 1.1)
        @NLconstraint(m, (x[1] - x[2])^2 >= 0.1)
        @NLconstraint(m, (x[2] - x[3])^2 >= 0.1)
        @NLconstraint(m, (x[1] - x[3])^2 >= 0.1)
        @NLconstraint(m, (x[1] - x[4])^2 >= 0.1)
        @NLconstraint(m, (x[2] - x[4])^2 >= 0.1)
        @NLconstraint(m, (x[3] - x[4])^2 >= 0.1)
        optimize!(m)
        nbranches = Juniper.getnbranches(unsafe_backend(m).inner)
        set_optimizer(
            m,
            optimizer_with_attributes(
                Juniper.Optimizer,
                juniper_all_solutions_p...,
            ),
        )
        JuMP.set_start_value.(x, zeros(4))
        optimize!(m)
        @test Juniper.getnbranches(unsafe_backend(m).inner) == nbranches
    end

    @testset "bruteforce PseudoCost" begin
        println("==================================")
        println("Bruteforce PseudoCost")
        println("==================================")
        juniper_all_solutions = DefaultTestSolver(
            branch_strategy = :PseudoCost,
            all_solutions = true,
            list_of_solutions = true,
            strong_restart = true,
            processors = 3,
        )
        m = Model(
            optimizer_with_attributes(
                Juniper.Optimizer,
                juniper_all_solutions...,
            ),
        )
        @variable(m, 1 <= x[1:4] <= 5, Int)
        @objective(m, Min, x[1])
        @constraint(m, x[1] >= 0.9)
        @constraint(m, x[1] <= 1.1)
        @NLconstraint(m, (x[1] - x[2])^2 >= 0.1)
        @NLconstraint(m, (x[2] - x[3])^2 >= 0.1)
        @NLconstraint(m, (x[1] - x[3])^2 >= 0.1)
        @NLconstraint(m, (x[1] - x[4])^2 >= 0.1)
        @NLconstraint(m, (x[2] - x[4])^2 >= 0.1)
        @NLconstraint(m, (x[3] - x[4])^2 >= 0.1)
        optimize!(m)
        status = termination_status(m)
        println("Status: ", status)
        list_of_solutions = Juniper.getsolutions(unsafe_backend(m).inner)
        @test length(unique(list_of_solutions)) == result_count(m)
        @test status == MOI.LOCALLY_SOLVED
        @test result_count(m) == 24
    end

    @testset "infeasible cos" begin
        println("==================================")
        println("Infeasible cos")
        println("==================================")
        m = Model(
            optimizer_with_attributes(
                Juniper.Optimizer,
                DefaultTestSolver(
                    branch_strategy = :StrongPseudoCost,
                    processors = 2,
                )...,
            ),
        )
        @variable(m, 1 <= x <= 5, Int)
        @variable(m, -2 <= y <= 2, Int)
        @objective(m, Min, -x - y)
        @NLconstraint(m, y == 2 * cos(2 * x))
        optimize!(m)
        status = termination_status(m)
        println("Status: ", status)
        @test status == MOI.LOCALLY_INFEASIBLE
    end

    @testset "infeasible relaxation" begin
        println("==================================")
        println("Infeasible relaxation")
        println("==================================")
        m = Model(
            optimizer_with_attributes(
                Juniper.Optimizer,
                DefaultTestSolver(
                    branch_strategy = :StrongPseudoCost,
                    processors = 2,
                )...,
            ),
        )
        @variable(m, 0 <= x[1:10] <= 2, Int)
        @objective(m, Min, sum(x))
        @constraint(m, sum(x[1:5]) <= 20)
        @NLconstraint(m, x[1] * x[2] * x[3] >= 10)
        optimize!(m)
        status = termination_status(m)
        println("Status: ", status)
        @test status == MOI.LOCALLY_INFEASIBLE
    end

    @testset "infeasible integer" begin
        println("==================================")
        println("Infeasible integer")
        println("==================================")
        m = Model(
            optimizer_with_attributes(
                Juniper.Optimizer,
                DefaultTestSolver(
                    branch_strategy = :StrongPseudoCost,
                    processors = 2,
                )...,
            ),
        )
        @variable(m, 0 <= x[1:10] <= 2, Int)
        @objective(m, Min, sum(x))
        @constraint(m, sum(x[1:5]) <= 20)
        @NLconstraint(m, x[1] * x[2] * x[3] >= 7)
        @NLconstraint(m, x[1] * x[2] * x[3] <= 7.5)
        optimize!(m)
        status = termination_status(m)
        println("Status: ", status)
        @test status == MOI.LOCALLY_INFEASIBLE
    end

    @testset "infeasible in strong" begin
        println("==================================")
        println("Infeasible in strong")
        println("==================================")
        m = Model(
            optimizer_with_attributes(
                Juniper.Optimizer,
                DefaultTestSolver(
                    branch_strategy = :StrongPseudoCost,
                    processors = 2,
                )...,
            ),
        )
        @variable(m, 0 <= x[1:5] <= 2, Int)
        @objective(m, Min, sum(x))
        @NLconstraint(m, x[3]^2 <= 2)
        @NLconstraint(m, x[3]^2 >= 1.2)
        optimize!(m)
        status = termination_status(m)
        println("Status: ", status)
        @test status == MOI.LOCALLY_INFEASIBLE
    end
end

include("POD_experiment/nous1.jl")
include("basic/gamsworld.jl")

@testset "POD instances" begin
    @testset "full strong branching" begin
        println("==================================")
        println("full strong branching")
        println("==================================")
        m = batch_problem()
        JuMP.set_optimizer(
            m,
            optimizer_with_attributes(
                Juniper.Optimizer,
                DefaultTestSolver(
                    branch_strategy = :StrongPseudoCost,
                    strong_branching_nsteps = 100,
                    strong_branching_perc = 100,
                    strong_restart = false,
                    debug = true,
                )...,
            ),
        )
        optimize!(m)
        @test termination_status(m) == MOI.LOCALLY_SOLVED
        juniper_val = JuMP.objective_value(m)
        println("obj: ", juniper_val)
        @test isapprox(
            juniper_val,
            285506.5082,
            atol = opt_atol,
            rtol = opt_rtol,
        )
        inner_model = unsafe_backend(m).inner
        counter_test(inner_model.debugDict, Juniper.getnbranches(inner_model))
    end

    @testset "break strong branching time limit" begin
        println("==================================")
        println("break strong branching time limit")
        println("==================================")
        m = batch_problem()
        set_optimizer(
            m,
            optimizer_with_attributes(
                Juniper.Optimizer,
                DefaultTestSolver(
                    branch_strategy = :StrongPseudoCost,
                    strong_branching_time_limit = 0.01,
                    time_limit = 4,
                    strong_restart = false,
                )...,
            ),
        )
        optimize!(m)
        status = termination_status(m)
        @test status == MOI.LOCALLY_SOLVED || status == MOI.TIME_LIMIT
        juniper_val = JuMP.objective_value(m)
        best_bound_val = JuMP.objective_bound(m)
        println("Best bound: ", best_bound_val)
        println("Objective: ", juniper_val)
        # minimization problem
        @test best_bound_val <= juniper_val || isnan(juniper_val)
    end

    @testset "nous1 restart" begin
        println("==================================")
        println("nous1 restart")
        println("==================================")
        m = get_nous1()
        set_optimizer(
            m,
            optimizer_with_attributes(
                Juniper.Optimizer,
                DefaultTestSolver(
                    branch_strategy = :StrongPseudoCost,
                    strong_restart = true,
                    mip_solver = optimizer_with_attributes(
                        HiGHS.Optimizer,
                        "output_flag" => false,
                    ),
                )...,
            ),
        )
        optimize!(m)
        if Sys.iswindows()
            # Spurious failure in CI: LOCALLY_INFEASIBLE
        else
            @test termination_status(m) == MOI.LOCALLY_SOLVED
        end
    end

    @testset "nous1 no restart" begin
        println("==================================")
        println("nous1 no restart")
        println("==================================")
        m = get_nous1()
        set_optimizer(
            m,
            optimizer_with_attributes(
                Juniper.Optimizer,
                DefaultTestSolver(
                    branch_strategy = :StrongPseudoCost,
                    strong_restart = false,
                    mip_solver = optimizer_with_attributes(
                        HiGHS.Optimizer,
                        "output_flag" => false,
                    ),
                )...,
            ),
        )
        optimize!(m)
        if Sys.iswindows()
            # Spurious failure in CI: LOCALLY_INFEASIBLE
        else
            @test termination_status(m) == MOI.LOCALLY_SOLVED
        end
    end

    @testset "reliability" begin
        println("==================================")
        println("reliability")
        println("==================================")
        m = batch_problem()
        set_optimizer(
            m,
            optimizer_with_attributes(
                Juniper.Optimizer,
                DefaultTestSolver(
                    branch_strategy = :Reliability,
                    reliability_branching_perc = 50,
                    reliability_branching_threshold = 5,
                    strong_restart = true,
                    debug = true,
                )...,
            ),
        )
        optimize!(m)
        @test termination_status(m) == MOI.LOCALLY_SOLVED
        juniper_val = JuMP.objective_value(m)
        best_bound_val = JuMP.objective_bound(m)
        gap_val = relative_gap(m)
        println("Solution by Juniper")
        println("obj: ", juniper_val)
        println("best_bound_val: ", best_bound_val)
        println("gap_val: ", gap_val)
        @test isapprox(juniper_val, 285506.5082, atol = 1e0)
        @test isapprox(best_bound_val, 285506.5082, atol = 1e0)
        @test isapprox(gap_val, 0, atol = 1e-2)
        inner_model = unsafe_backend(m).inner
        counter_test(inner_model.debugDict, Juniper.getnbranches(inner_model))
    end
end

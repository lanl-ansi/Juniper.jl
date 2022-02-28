include("POD_experiment/FLay02H.jl")

@testset "User limit testing" begin
    @testset "Knapsack 50% limit" begin
        println("==================================")
        println("KNAPSACK 50%")
        println("==================================")
        m = Model(
            optimizer_with_attributes(
                Juniper.Optimizer,
                DefaultTestSolver(;
                    traverse_strategy = :DBFS,
                    branch_strategy = :MostInfeasible,
                    mip_gap = 0.5,
                )...,
            ),
        )
        v = [10, 20, 12, 23, 42]
        w = [12, 45, 12, 22, 21]
        @variable(m, x[1:5], Bin, start = 0.5)
        @objective(m, Max, v' * x)
        @NLconstraint(m, sum(w[i] * x[i]^2 for i in 1:5) <= 45)
        optimize!(m)
        objval = JuMP.objective_value(m)
        println("Obj: ", objval)
        @test termination_status(m) == MOI.OBJECTIVE_LIMIT
        @test JuMP.objective_bound(m) >= objval
        @test 0.1 <= relative_gap(m) <= 0.5
    end

    @testset "FLay02H time limit" begin
        println("==================================")
        println("FLay02H time limit")
        println("==================================")
        m = get_FLay02H()
        set_optimizer(
            m,
            optimizer_with_attributes(
                Juniper.Optimizer,
                DefaultTestSolver(
                    branch_strategy = :StrongPseudoCost,
                    time_limit = 5,
                    incumbent_constr = true,
                )...,
            ),
        )
        optimize!(m)
        status = termination_status(m)
        @test status == MOI.LOCALLY_SOLVED || status == MOI.TIME_LIMIT
        @test solve_time(m) <= 15 # it might be a bit higher than 5s
    end
end

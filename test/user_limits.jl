include("POD_experiment/blend029.jl")

@testset "User limit testing" begin

@testset "Knapsack 50% limit" begin
    println("==================================")
    println("KNAPSACK 50%")
    println("==================================")

    m = Model(with_optimizer(
        Juniper.Optimizer,
        DefaultTestSolver(;traverse_strategy=:DBFS,branch_strategy=:MostInfeasible,mip_gap=0.5))
    )

    v = [10,20,12,23,42]
    w = [12,45,12,22,21]
    @variable(m, x[1:5], Bin)

    @objective(m, Max, dot(v,x))

    @NLconstraint(m, sum(w[i]*x[i]^2 for i=1:5) <= 45)   

    status = solve(m)
    objval = JuMP.objective_value(m)
    println("Obj: ", objval)
    best_bound_val = JuMP.objective_bound(m)
    gap_val = getobjgap(m)

    @test status == MOI.OBJECTIVE_LIMIT

    @test best_bound_val >= objval
    @test 0.1 <= gap_val <= 0.5
end

@testset "blend029 1s limit" begin
    println("==================================")
    println("blend029 1s limit")
    println("==================================")

    m,objval = get_blend029()

    set_optimizer(m, with_optimizer(
        Juniper.Optimizer,
        DefaultTestSolver(
            branch_strategy=:StrongPseudoCost, 
            time_limit = 1, # second
            incumbent_constr = true
        )
    ))

    status = solve(m)

    @test status == MOI.TIME_LIMIT

    juniper_val = JuMP.objective_value(m)
    best_bound_val = JuMP.objective_bound(m)
    gap_val = getobjgap(m)

    println("Solution by Juniper")
    println("obj: ", juniper_val)
    println("best_bound_val: ", best_bound_val)
    println("gap_val: ", gap_val)

    @test best_bound_val >= objval
    @test getsolvetime(m) <= 4 # it might be a bit higher than 1s
end

end
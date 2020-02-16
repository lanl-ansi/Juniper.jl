include("POD_experiment/FLay02H.jl")

@testset "User limit testing" begin

@testset "Knapsack 50% limit" begin
    println("==================================")
    println("KNAPSACK 50%")
    println("==================================")

    m = Model(optimizer_with_attributes(
        Juniper.Optimizer,
        DefaultTestSolver(;traverse_strategy=:DBFS,branch_strategy=:MostInfeasible,mip_gap=0.5)...)
    )

    v = [10,20,12,23,42]
    w = [12,45,12,22,21]
    @variable(m, x[1:5], Bin, start=0.5)

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

@testset "FLay02H time limit" begin
    println("==================================")
    println("FLay02H time limit")
    println("==================================")

    m = get_FLay02H()

    set_optimizer(m, optimizer_with_attributes(
        Juniper.Optimizer,
        DefaultTestSolver(
            branch_strategy=:StrongPseudoCost,
            time_limit = 5,
            incumbent_constr = true)...
    ))

    status = solve(m)

    @test status == MOI.LOCALLY_SOLVED || status == MOI.TIME_LIMIT
    @test getsolvetime(m) <= 15 # it might be a bit higher than 5s
end

end
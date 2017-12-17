include("POD_experiment/blend029.jl")

@testset "User limit testing" begin

@testset "Knapsack 50% limit" begin
    println("==================================")
    println("KNAPSACK 50%")
    println("==================================")

    m = Model(solver=DefaultTestSolver(;traverse_strategy=:DBFS,branch_strategy=:MostInfeasible,mip_gap=0.5))

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

    @test status == :UserLimit

    @test best_bound_val >= objval
    @test 0.1 <= gap_val <= 0.5
end

@testset "blend029 10s limit" begin
    println("==================================")
    println("blend029 10s limit")
    println("==================================")

    m,objval = get_blend029()

    setsolver(m, DefaultTestSolver(
            branch_strategy=:PseudoCost, 
            time_limit = 10, # seconds
            incumbent_constr = true
    ))
    status = solve(m)

    @test status == :UserLimit

    juniper_val = getobjectivevalue(m)
    best_bound_val = getobjbound(m)
    gap_val = getobjgap(m)

    println("Solution by Juniper")
    println("obj: ", juniper_val)
    println("best_bound_val: ", best_bound_val)
    println("gap_val: ", gap_val)

    @test best_bound_val >= objval
    @test getsolvetime(m) <= 15 # it might be a bit higher than 10s
end

@testset "case 5 socwr solution limit" begin
    println("==============================================")
    println("SOCWRPowerModel case5.m SolutionLimit=1")
    println("==============================================")

    pm = build_generic_model("data/pglib_opf_case5_pjm.m", SOCWRPowerModel, PowerModels.post_ots)
    m = pm.model
    @variable(m, 0 <= aeiou <= 1)
    @NLconstraint(m, aeiou^2== 1)

    solver = DefaultTestSolver(
        branch_strategy=:StrongPseudoCost, #
        strong_restart = false,
        solution_limit = 1
    )
    setsolver(m, solver)
    status = solve(m)

    @test status == :UserLimit
    nsolutions = Juniper.getnsolutions(internalmodel(m))

    juniper_val = getobjectivevalue(m)

    println("Solution by Juniper")
    println("obj: ", juniper_val)

    @test nsolutions == 1 || nsolutions == 2 # if last branch has two integral solutions
end

@testset "case 5 socwr best_obj_stop" begin
    println("==============================================")
    println("SOCWRPowerModel case5.m best_obj_stop=15000")
    println("==============================================")

    pm = build_generic_model("data/pglib_opf_case5_pjm.m", SOCWRPowerModel, PowerModels.post_ots)
    m = pm.model
    @variable(m, 0 <= aeiou <= 1)
    @NLconstraint(m, aeiou^2== 1)

    best_obj_stop = 15000
    solver = DefaultTestSolver(
        branch_strategy=:StrongPseudoCost, #
        strong_restart = false,
        best_obj_stop = best_obj_stop
    )
    setsolver(m, solver)
    status = solve(m)

    @test status == :UserLimit
    nsolutions = Juniper.getnsolutions(internalmodel(m))

    juniper_val = getobjectivevalue(m)

    println("Solution by Juniper")
    println("obj: ", juniper_val)

    @test juniper_val <= best_obj_stop
end

end
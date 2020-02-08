@testset "Power Models SOCWR tests" begin

@testset "case 5" begin
    println("==================================")
    println("SOCWRPowerModel case5.m")
    println("==================================")

    pm = build_generic_model("data/pglib_opf_case5_pjm.m", SOCWRPowerModel, PowerModels.post_ots)
    m = pm.model
    @variable(m, 0 <= aeiou <= 1)
    @NLconstraint(m, aeiou^2== 1)

    setsolver(m, juniper_strong_no_restart)
    status = solve(m)

    @test status == :Optimal || status == :LocalOptimal

    juniper_val = getobjectivevalue(m)

    println("Solution by Juniper")
    println("obj: ", juniper_val)

    @test isapprox(juniper_val, 14999.7, atol=1e0)
end

@testset "case 5 w inc constr" begin
    println("==============================================")
    println("SOCWRPowerModel case5.m no incumbent constr")
    println("==============================================")

    pm = build_generic_model("data/pglib_opf_case5_pjm.m", SOCWRPowerModel, PowerModels.post_ots)
    m = pm.model
    @variable(m, 0 <= aeiou <= 1)
    @NLconstraint(m, aeiou^2== 1)

    solver = DefaultTestSolver(
        branch_strategy=:StrongPseudoCost,
        strong_restart = false,
        incumbent_constr = true,
        traverse_strategy = :DFS
    )
    setsolver(m, solver...)
    status = solve(m)

    @test status == :Optimal || status == :LocalOptimal

    juniper_val = getobjectivevalue(m)

    println("Solution by Juniper")
    println("obj: ", juniper_val)

    @test m.internalModel.options.strong_restart == false
    @test m.internalModel.options.incumbent_constr == true
    @test isapprox(juniper_val, 14999.7, atol=1e0)
end

@testset "case 5 epsilon obj constr" begin
    println("==============================================")
    println("SOCWRPowerModel case5.m epsilon constr")
    println("==============================================")

    pm = build_generic_model("data/pglib_opf_case5_pjm.m", SOCWRPowerModel, PowerModels.post_ots)
    m = pm.model
    @variable(m, 0 <= aeiou <= 1)
    @NLconstraint(m, aeiou^2== 1)

    solver = DefaultTestSolver(
        branch_strategy=:StrongPseudoCost,
        strong_branching_perc = 0, # this defaults to 2 variables
        strong_restart = false,
        obj_epsilon = 0.5,
        traverse_strategy = :DFS
    )
    setsolver(m, solver...)
    status = solve(m)

    @test status == :Optimal || status == :LocalOptimal

    juniper_val = getobjectivevalue(m)

    println("Solution by Juniper")
    println("obj: ", juniper_val)

    @test isapprox(juniper_val, 14999.7, atol=1e0)
end

end
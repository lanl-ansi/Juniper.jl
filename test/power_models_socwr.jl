@testset "Power Models SOCWR tests" begin

@testset "case 5" begin
    println("==================================")
    println("SOCWRPowerModel case5.m")
    println("==================================")

    pm = build_generic_model("data/pglib_opf_case5_pjm.m", SOCWRPowerModel, PowerModels.post_ots)
    m = pm.model
    @variable(m, 0 <= aeiou <= 1)
    @NLconstraint(m, aeiou^2== 1)

    setsolver(m, minlpbnb_strong_no_restart)
    status = solve(m)

    @test status == :Optimal

    minlpbnb_val = getobjectivevalue(m)

    println("Solution by MINLPBnb")
    println("obj: ", minlpbnb_val)

    @test isapprox(minlpbnb_val, 14999.7, atol=1e0)
end

@testset "case 5 wo inc constr" begin
    println("==============================================")
    println("SOCWRPowerModel case5.m no incumbent constr")
    println("==============================================")

    pm = build_generic_model("data/pglib_opf_case5_pjm.m", SOCWRPowerModel, PowerModels.post_ots)
    m = pm.model
    @variable(m, 0 <= aeiou <= 1)
    @NLconstraint(m, aeiou^2== 1)

    solver = DefaultTestSolver(
        branch_strategy=:StrongPseudoCost,
        strong_branching_nvars = 5,
        strong_restart = false,
        incumbent_constr = false,
        traverse_strategy = :DFS
    )
    setsolver(m, solver)
    status = solve(m)

    @test status == :Optimal

    minlpbnb_val = getobjectivevalue(m)

    println("Solution by MINLPBnb")
    println("obj: ", minlpbnb_val)

    @test isapprox(minlpbnb_val, 14999.7, atol=1e0)
end

end
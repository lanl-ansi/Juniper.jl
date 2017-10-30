@testset "Power Models ACP tests" begin

@testset "case 5" begin
    println("==================================")
    println("ACPPowerModel case5.m")
    println("==================================")

    result = run_ots("data/pglib_opf_case5_pjm.m", ACPPowerModel, minlpbnb_strong_no_restart)

    status = result["status"]
    @test status == :Optimal

    minlpbnb_val = result["objective"]

    println("Solution by MINLPBnb")
    println("obj: ", minlpbnb_val)

    @test isapprox(minlpbnb_val, 15174.0, atol=1e0)
end

@testset "case14" begin
    println("==================================")
    println("ACPPowerModel case14.m")
    println("==================================")

    result = run_ots("data/pglib_opf_case14_ieee.m", ACPPowerModel, minlpbnb_strong_no_restart)

    status = result["status"]
    @test status == :Optimal

    minlpbnb_val = result["objective"]

    println("")
    println("Solution by MINLPBnb")
    println("obj: ", minlpbnb_val)

    @test isapprox(minlpbnb_val, 6291.28, atol=1e0)
end

end
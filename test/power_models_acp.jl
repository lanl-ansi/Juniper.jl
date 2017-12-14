@testset "Power Models ACP tests" begin

@testset "case14" begin
    println("==================================")
    println("ACPPowerModel case14.m")
    println("==================================")

    result = run_ots("data/pglib_opf_case14_ieee.m", ACPPowerModel, juniper_strong_no_restart)

    status = result["status"]
    @test status == :Optimal || status == :LocalOptimal

    juniper_val = result["objective"]

    println("")
    println("Solution by Juniper")
    println("obj: ", juniper_val)

    @test isapprox(juniper_val, 6291.28, atol=1e0)
end

end
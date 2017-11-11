include("POD_experiment/blend029.jl")

@testset "POD instances" begin

@testset "blend029 full strong branching" begin
    println("==================================")
    println("blend029 full strong branching")
    println("==================================")

    m,objval = get_blend029()

    setsolver(m, DefaultTestSolver(
            branch_strategy=:StrongPseudoCost,
            strong_branching_perc = 100,
            strong_branching_nsteps = 100,
            strong_restart = true
    ))
    status = solve(m)

    @test status == :Optimal

    minlpbnb_val = getobjectivevalue(m)
    best_bound_val = getobjbound(m)
    gap_val = getobjgap(m)

    println("Solution by MINLPBnb")
    println("obj: ", minlpbnb_val)
    println("best_bound_val: ", best_bound_val)
    println("gap_val: ", gap_val)

    @test isapprox(minlpbnb_val, objval, atol=1e0)
    @test isapprox(best_bound_val, objval, atol=1e0)
    @test isapprox(gap_val, 0, atol=1e-2)
end


end
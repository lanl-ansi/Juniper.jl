include("POD_experiment/blend029.jl")

@testset "POD instances" begin

@testset "blend029" begin
    println("==================================")
    println("blend029")
    println("==================================")

    m,objval = get_blend029()

    setsolver(m, MINLPBnBSolver(IpoptSolver(print_level=0);
            branch_strategy=:StrongPseudoCost,
            strong_branching_nvars = 15,
            strong_branching_nsteps = 5,
            strong_restart = true
    ))
    status = solve(m)

    @test status == :Optimal

    minlpbnb_val = getobjectivevalue(m)

    println("Solution by MINLPBnb")
    println("obj: ", minlpbnb_val)

    @test isapprox(minlpbnb_val, objval, atol=1e0)
end

end
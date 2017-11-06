include("basic/gamsworld.jl")
include("POD_experiment/blend029.jl")

@testset "parallel tests" begin

@testset "Batch.mod no restart parallel" begin
    println("==================================")
    println("BATCH.MOD NO RESTART")
    println("==================================")

    m = batch_problem()

    minlpbnb = MINLPBnBSolver(IpoptSolver(print_level=0);
        branch_strategy=:StrongPseudoCost,
        strong_restart = false,
        processors = 4
    ) 

    setsolver(m, minlpbnb)
    status = solve(m)
    @test status == :Optimal

    minlpbnb_val = getobjectivevalue(m)
    minlpbnb_bb = getobjbound(m)

    println("Solution by MINLPBnb")
    println("obj: ", minlpbnb_val)
    println("bound: ", minlpbnb_bb)


    @test isapprox(minlpbnb_val, 285506.5082, atol=opt_atol, rtol=opt_rtol)
end

@testset "blend029" begin
    println("==================================")
    println("blend029")
    println("==================================")

    m,objval = get_blend029()

    setsolver(m, MINLPBnBSolver(IpoptSolver(print_level=0);
            branch_strategy=:StrongPseudoCost,
            strong_branching_nvars = 15,
            strong_branching_nsteps = 5,
            strong_restart = true,
            processors = 4
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
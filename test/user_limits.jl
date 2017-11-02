include("POD_experiment/blend029.jl")

@testset "User limit testing" begin

@testset "blend029 10s limit" begin
    println("==================================")
    println("blend029 10s limit")
    println("==================================")

    m,objval = get_blend029()

    setsolver(m, MINLPBnBSolver(IpoptSolver(print_level=0);
            branch_strategy=:StrongPseudoCost,
            strong_branching_nvars = 15,
            strong_branching_nsteps = 5,
            strong_restart = true,
            time_limit = 10 # seconds
    ))
    status = solve(m)

    @test status == :UserLimit

    minlpbnb_val = getobjectivevalue(m)
    best_bound_val = getobjbound(m)
    gap_val = getobjgap(m)

    println("Solution by MINLPBnb")
    println("obj: ", minlpbnb_val)
    println("best_bound_val: ", best_bound_val)
    println("gap_val: ", gap_val)

    @test best_bound_val >= objval
end

@testset "blend029 5% limit" begin
    println("==================================")
    println("blend029 5%")
    println("==================================")

    m,objval = get_blend029()

    setsolver(m, MINLPBnBSolver(IpoptSolver(print_level=0);
            branch_strategy=:StrongPseudoCost,
            strong_branching_nvars = 15,
            strong_branching_nsteps = 5,
            strong_restart = true,
            mip_gap = 5 # %
    ))
    status = solve(m)

    @test status == :UserLimit

    minlpbnb_val = getobjectivevalue(m)
    best_bound_val = getobjbound(m)
    gap_val = getobjgap(m)

    println("Solution by MINLPBnb")
    println("obj: ", minlpbnb_val)
    println("best_bound_val: ", best_bound_val)
    println("gap_val: ", gap_val)

    @test best_bound_val >= objval
    @test gap_val <= 0.05
end

@testset "case 5 socwr solution limit" begin
    println("==============================================")
    println("SOCWRPowerModel case5.m SolutionLimit=1")
    println("==============================================")

    pm = build_generic_model("data/pglib_opf_case5_pjm.m", SOCWRPowerModel, PowerModels.post_ots)
    m = pm.model
    @variable(m, 0 <= aeiou <= 1)
    @NLconstraint(m, aeiou^2== 1)

    solver = MINLPBnBSolver(IpoptSolver(print_level=0);
        branch_strategy=:StrongPseudoCost,
        strong_branching_nvars = 5,
        strong_restart = false,
        incumbent_constr = false,
        solution_limit = 1
    )
    setsolver(m, solver)
    status = solve(m)

    @test status == :UserLimit
    nsolutions = MINLPBnB.getnsolutions(internalmodel(m))

    minlpbnb_val = getobjectivevalue(m)

    println("Solution by MINLPBnb")
    println("obj: ", minlpbnb_val)

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
    solver = MINLPBnBSolver(IpoptSolver(print_level=0);
        branch_strategy=:StrongPseudoCost,
        strong_branching_nvars = 5,
        strong_restart = false,
        incumbent_constr = false,
        best_obj_stop = best_obj_stop
    )
    setsolver(m, solver)
    status = solve(m)

    @test status == :UserLimit
    nsolutions = MINLPBnB.getnsolutions(internalmodel(m))

    minlpbnb_val = getobjectivevalue(m)

    println("Solution by MINLPBnb")
    println("obj: ", minlpbnb_val)

    @test minlpbnb_val <= best_obj_stop
end

end
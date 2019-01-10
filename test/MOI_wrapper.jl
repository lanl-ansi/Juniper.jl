import MathOptInterface
const MOI = MathOptInterface
const MOIT = MOI.Test

const optimizer = Juniper.Optimizer(
    solver = DefaultTestSolver(
        branch_strategy=:StrongPseudoCost,
        strong_branching_perc = 25,
        strong_branching_nsteps = 2,
        strong_restart = true
    )
)

const config = MOIT.TestConfig(
    atol = 1e-4, rtol = 1e-4, optimal_status = MOI.LOCALLY_SOLVED
)

@testset "MOI NLP tests" begin
    MOIT.nlptest(optimizer, config)
end

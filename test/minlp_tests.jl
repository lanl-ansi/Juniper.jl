using MINLPTests

const OPTIMIZER = MINLPTests.JuMP.optimizer_with_attributes(
    Juniper.Optimizer, DefaultTestSolver()...
)

@testset "MINLPTests" begin
    MINLPTests.test_nlp_mi(OPTIMIZER, exclude = [
        "005_011",  # Uses the function `\`
        "003_013",  # Bug in Juniper - handling of expression graph?
        "006_010"   # Bug in Juniper - handling of user-defined functions.
    ], objective_tol = 1e-5, primal_tol = 1e-5, dual_tol = NaN)
end
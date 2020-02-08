using MINLPTests, JuMP, Ipopt, Juniper, Test

const OPTIMIZER = MINLPTests.JuMP.optimizer_with_attributes(
    Juniper.Optimizer, nl_solver=optimizer_with_attributes(Ipopt.Optimizer, print_level=0), atol=1e-7
)

@testset "MINLPTests" begin
    ###
    ### src/nlp-mi tests.
    ###
    
    MINLPTests.test_nlp_mi(OPTIMIZER)
end
using MINLPTests, JuMP, Ipopt, Juniper, Test

const OPTIMIZER = MINLPTests.JuMP.with_optimizer(
    Juniper.Optimizer, nl_solver=with_optimizer(Ipopt.Optimizer, print_level=0), atol=1e-7
)

@testset "MINLPTests" begin
    ###
    ### src/nlp-mi tests.
    ###
    
    MINLPTests.test_nlp_mi(OPTIMIZER)
end
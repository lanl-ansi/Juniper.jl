module TestMOIWrapper

using Juniper
using MathOptInterface
using Test
import Ipopt

const MOI = MathOptInterface

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$(name)", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
    return
end

function test_SolverName()
    @test MOI.get(Juniper.Optimizer(), MOI.SolverName()) == "Juniper"
    return
end

function test_runtests()
    config = MOI.Test.Config(
        atol = 1e-4,
        rtol = 1e-4,
        optimal_status = MOI.LOCALLY_SOLVED,
        exclude = Any[
            MOI.ConstraintDual,
            MOI.DualObjectiveValue,
            MOI.ConstraintBasisStatus,
            MOI.VariableBasisStatus,
            MOI.NLPBlockDual,
        ],
    )
    optimizer = Juniper.Optimizer(
        Pair{String, Any}[
            "log_levels" => [:Table],
            "nl_solver" => MOI.OptimizerWithAttributes(
                Ipopt.Optimizer,
                "print_level" => 0,
                "sb" => "yes",
            ),
        ],
    )
    caching_optimizer = MOI.Utilities.CachingOptimizer(
        MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
        optimizer,
    )
    MOI.Test.runtests(
        caching_optimizer,
        config;
        exclude = String[
            # ====================== Unexplained failures ======================
            # LOCALLY_SOLVED instead of INFEASIBLE
            "test_conic_linear_INFEASIBLE_2",
            # NORM_LIMIT instead of LOCALLY_SOLVED
            "test_conic_linear_VectorOfVariables",
            # INVALID_MODEL instead of LOCALLY_SOLVED
            "test_linear_VectorAffineFunction_empty_row",
            # Attributes are not checked properly on copy_to
            "test_model_copy_to_UnsupportedAttribute",
            "test_model_copy_to_UnsupportedConstraint",
            # supports_constraint not checking for element types
            "test_model_supports_constraint_ScalarAffineFunction_EqualTo",
            "test_model_supports_constraint_VariableIndex_EqualTo",
            "test_model_supports_constraint_VectorOfVariables_Nonnegatives",
            # Result index not being checked properly
            "test_solve_result_index",
            # ======================= Explained failures =======================
            # LOCALLY_INFEASIBLE instead of INFEASIBLE
            "test_conic_NormInfinityCone_INFEASIBLE",
            "test_conic_NormOneCone_INFEASIBLE",
            "test_conic_linear_INFEASIBLE",
            "test_constraint_ZeroOne_bounds_3",
            "test_linear_INFEASIBLE",
            "test_linear_INFEASIBLE_2",
            # NORM_LIMIT instead of DUAL_INFEASIBLE
            "test_solve_TerminationStatus_DUAL_INFEASIBLE",
            "test_linear_DUAL_INFEASIBLE",
            "test_linear_DUAL_INFEASIBLE_2",
        ]
    )
    return
end

end

TestMOIWrapper.runtests()

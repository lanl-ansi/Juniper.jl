using MathOptInterface
const MOI = MathOptInterface
const MOIT = MOI.Test
const MOIU = MOI.Utilities
const MOIB = MOI.Bridges

juniper = DefaultTestSolver(log_levels=[:Table])

const optimizer = Juniper.Optimizer(juniper)

@testset "SolverName" begin
    @test MOI.get(optimizer, MOI.SolverName()) == "Juniper"
end

const config = MOIT.Config(atol=1e-4, rtol=1e-4, optimal_status=MOI.LOCALLY_SOLVED, exclude = Any[
    MOI.ConstraintDual,
    MOI.VariableName,
    MOI.ConstraintName,
    MOI.delete,
])

@testset "Unit" begin
    bridged = MOIB.full_bridge_optimizer(optimizer, Float64)
    # A number of test cases are excluded because loadfromstring! works only
    # if the solver supports variable and constraint names.
    exclude = [
                "delete_nonnegative_variables", # delete and ConstraintFunction not supported
                "delete_soc_variables", # Soc and delete not supported
                "delete_variable", # Deleting not supported.
                "delete_variables", # Deleting not supported.
                "getvariable", # Variable names not supported.
                "solve_zero_one_with_bounds_1", # Variable names not supported.
                "solve_zero_one_with_bounds_2", # Variable names not supported.
                "solve_zero_one_with_bounds_3", # Variable names not supported.
                "getconstraint", # Constraint names not suported.
                "variablenames", # Variable names not supported.
                "solve_with_upperbound", # loadfromstring!
                "solve_with_lowerbound", # loadfromstring!
                "solve_integer_edge_cases", # loadfromstring!
                "solve_affine_lessthan", # loadfromstring!
                "solve_affine_greaterthan", # loadfromstring!
                "solve_affine_equalto", # loadfromstring!
                "solve_affine_interval", # loadfromstring!
                "get_objective_function", # Function getters not supported.
                "solve_constant_obj",  # loadfromstring!
                "solve_blank_obj", # loadfromstring!
                "solve_singlevariable_obj", # loadfromstring!
                "solve_objbound_edge_cases", # ObjectiveBound not supported.
                "solve_affine_deletion_edge_cases", # Deleting not supported.
                "solve_unbounded_model", # `NORM_LIMIT`
                "solve_duplicate_terms_obj", # duplicate terms objective
                "solve_duplicate_terms_vector_affine", 
                "solve_duplicate_terms_scalar_affine", 
                "solve_result_index", # no support for `MOI.ConstraintPrimal` atm
                "update_dimension_nonnegative_variables", # currently no support for ConstraintFunction 
               ]
    MOIT.runtests(bridged, config, exclude)
end

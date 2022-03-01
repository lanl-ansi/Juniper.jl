module TestConicModels

using Test
using JuMP

import Juniper
import SCS

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

function test_simple_conic_model()
    model = Model(
        optimizer_with_attributes(
            Juniper.Optimizer,
            "nl_solver" =>
                optimizer_with_attributes(SCS.Optimizer, MOI.Silent() => true),
            "atol" => 1e-4,
        ),
    )
    @variable(model, 0 <= x[1:2] <= 10, Int)
    @objective(model, Max, 3x[1] + 5x[2])
    @constraint(model, [10, x[1], x[2]] in SecondOrderCone())
    optimize!(model)
    @test termination_status(model) == LOCALLY_SOLVED
    @test primal_status(model) == FEASIBLE_POINT
    @test dual_status(model) == FEASIBLE_POINT
    @test isapprox(value.(x), [6, 8]; atol = 1e-4)
    @test sqrt(value(x[1])^2 + value(x[2])^2) <= 10 + 1e-4
    @test isapprox(objective_value(model), 58; atol = 1e-5)
    return
end

end

TestConicModels.runtests()

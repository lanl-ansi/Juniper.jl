
function DefaultTestSolver(;nl_solver=optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0, "sb" =>"yes"), solver_args...)
    solver_args_result = Vector{Pair{String, Any}}()
    push!(solver_args_result, "log_levels" => Symbol[])
    push!(solver_args_result, "nl_solver" => nl_solver)
    for v in solver_args
        push!(solver_args_result, string(v[1]) => v[2])
    end
    return solver_args_result
end
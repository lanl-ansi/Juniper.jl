export MINLPBnBSolver

"""
A solver for MINLP problems using a NLP solver and Branch and Bound
"""

type MINLPBnBSolverObj <: MathProgBase.AbstractMathProgSolver
    nl_solver   :: MathProgBase.AbstractMathProgSolver
    options     :: MINLPBnB.SolverOptions
end

function MINLPBnBSolver(nl_solver::MathProgBase.AbstractMathProgSolver;
        log_levels                  = [:Table,:Info],
        branch_strategy             = :MostInfeasible,
        strong_branching_nvars      = 5,
        strong_branching_nlevels    = 1
    )
    options_obj = MINLPBnB.SolverOptions(log_levels,
                                        branch_strategy,
                                        strong_branching_nvars,
                                        strong_branching_nlevels)
    return MINLPBnBSolverObj(nl_solver,options_obj)
end
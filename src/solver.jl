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
        strong_branching_nsteps     = 1,
        strong_restart              = true,
        incumbent_constr            = true
    )
    options_obj = MINLPBnB.SolverOptions(log_levels,
                                        branch_strategy,
                                        strong_branching_nvars,
                                        strong_branching_nsteps,
                                        strong_restart,
                                        incumbent_constr)
    return MINLPBnBSolverObj(nl_solver,options_obj)
end
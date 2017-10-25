export MINLPBnBSolver

"""
A solver for MINLP problems using a NLP solver and Branch and Bound
"""

type MINLPBnBSolverObj <: MathProgBase.AbstractMathProgSolver
    nl_solver   :: MathProgBase.AbstractMathProgSolver
    print_syms  :: Vector{Symbol}
end

function MINLPBnBSolver(nl_solver::MathProgBase.AbstractMathProgSolver;print_syms=[])
    return MINLPBnBSolverObj(nl_solver,print_syms)
end
export MINLPBnBSolver

"""
A solver for MINLP problems using a NLP solver and Branch and Bound
"""

type MINLPBnBSolverObj <: MathProgBase.AbstractMathProgSolver
    nl_solver   :: MathProgBase.AbstractMathProgSolver
end

function MINLPBnBSolver(nl_solver::MathProgBase.AbstractMathProgSolver)
    println("solver.jl => MINLPBnB")
    return MINLPBnBSolverObj(nl_solver)
end
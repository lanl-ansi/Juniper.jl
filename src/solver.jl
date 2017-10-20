export MINLPBnBSolver

"""
A solver for MIQCQP problems using a QCQP solver and Branch and Bound
"""
type MINLPBnBSolverObj <: MathProgBase.AbstractMathProgSolver
    nl_solver   :: MathProgBase.AbstractMathProgSolver
end

function MINLPBnBSolver(nl_solver::MathProgBase.AbstractMathProgSolver)
    println("solver.jl => MINLPBnB")
    return MINLPBnBSolverObj(nl_solver)
end

MathProgBase.LinearQuadraticModel(s::MINLPBnBSolverObj) = MathProgBase.NonlinearToLPQPBridge(MathProgBase.NonlinearModel(s))
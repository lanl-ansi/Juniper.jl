import MathOptInterface

const MOI  = MathOptInterface
const MOIU = MathOptInterface.Utilities


# JuniperProblem
mutable struct JuniperProblem
    nl_solver       :: MOI.AbstractOptimizer
       
    model           :: JuMP.Model
            
    status          :: Symbol
    objval          :: Float64
    best_bound      :: Float64
    
    x               :: Vector{JuMP.VariableRef}
    num_constr      :: Int64
    num_nl_constr   :: Int64
    num_l_constr    :: Int64
    num_var         :: Int64
    l_var           :: Vector{Float64}
    u_var           :: Vector{Float64}
    l_constr        :: Vector{Float64}
    u_constr        :: Vector{Float64}
    
    affs            :: Vector{Aff}
    
    disc2var_idx    :: Vector{Int64}
    var2disc_idx    :: Vector{Int64}
    
    var_type        :: Vector{Symbol}
    isconstrlinear  :: Vector{Bool}
    obj_sense       :: Symbol
    d               :: MathProgBase.AbstractNLPEvaluator
    num_disc_var :: Int64
    
    solution        :: Vector{Float64}
    
    soltime         :: Float64
    options         :: SolverOptions
    solutions       :: Vector{SolutionObj}
    nsolutions      :: Int64
    
    mip_solver      :: MathProgBase.AbstractMathProgSolver
    
    relaxation_time :: Float64
    start_time      :: Float64
    
    # Info
    nintvars        :: Int64
    nbinvars        :: Int64
    nnodes          :: Int64
    ncuts           :: Int64
    nbranches       :: Int64
    nlevels         :: Int64
    
    fpump_info      :: Dict{Symbol,Float64}
    
    # debug
    debugDict        :: Dict{Symbol,Any}
    
    JuniperProblem() = new()
end


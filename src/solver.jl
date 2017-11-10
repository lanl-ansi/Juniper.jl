export MINLPBnBSolver

"""
A solver for MINLP problems using a NLP solver and Branch and Bound
"""

type MINLPBnBSolverObj <: MathProgBase.AbstractMathProgSolver
    nl_solver   :: MathProgBase.AbstractMathProgSolver
    options     :: MINLPBnB.SolverOptions
end

function get_default_options()
    log_levels                  = [:Table,:Info]
    branch_strategy             = :StrongPseudoCost
    # Strong branching
    strong_branching_nvars      = 5
    strong_branching_nsteps     = 1
    strong_restart              = true
    # Obj cuts
    incumbent_constr            = true
    obj_epsilon                 = 0
    # :UserLimit
    time_limit                  = Inf  
    mip_gap                     = 1e-2 # in % used in bb_user_limits (#TODO)
    best_obj_stop               = NaN
    solution_limit              = 0
    all_solutions               = false
    list_of_solutions           = false
    # Parallel
    processors                  = 1
    # Traversing
    traverse_strategy           = :BFS
    return SolverOptions(log_levels,branch_strategy,strong_branching_nvars,strong_branching_nsteps,strong_restart,
        incumbent_constr,obj_epsilon,time_limit,mip_gap,best_obj_stop,solution_limit,all_solutions,
        list_of_solutions,processors,traverse_strategy)
end

function combine_options(options)
    options_dict = Dict{Symbol,Any}()
    for kv in options
        options_dict[kv[1]] = kv[2]
    end
    if get(options_dict,:log_levels,false) != false
        if length(options_dict[:log_levels]) == 0
            options_dict[:log_levels] = Symbol[]
        end
    end
    defaults = get_default_options()
    for fname in fieldnames(SolverOptions)
        if get(options_dict,fname,false) != false
            if fieldtype(SolverOptions,fname) != typeof(options_dict[fname])
                options_dict[fname] = convert(fieldtype(SolverOptions,fname),options_dict[fname])
            end
            setfield!(defaults,fname,options_dict[fname])
        end
    end
    return defaults
end

function MINLPBnBSolver(nl_solver::MathProgBase.AbstractMathProgSolver;options...)
    options_obj = combine_options(options)
    return MINLPBnBSolverObj(nl_solver,options_obj)
end
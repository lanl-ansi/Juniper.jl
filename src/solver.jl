export MINLPBnBSolver

"""
A solver for MINLP problems using a NLP solver and Branch and Bound
"""

type MINLPBnBSolverObj <: MathProgBase.AbstractMathProgSolver
    nl_solver   :: MathProgBase.AbstractMathProgSolver
    options     :: MINLPBnB.SolverOptions
end

function get_default_options()
    log_levels                      = [:Options,:Table,:Info]
    branch_strategy                 = :StrongPseudoCost
    gain_mu                         = 0.167 # http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.92.7117&rep=rep1&type=pdf
    # Strong branching
    strong_branching_perc           = 25
    strong_branching_nsteps         = 1
    strong_restart                  = true
    # Reliability branching 
    reliability_branching_threshold = 5 # reliability param
    reliability_branching_perc      = 25
    # Obj cuts
    incumbent_constr                = true
    obj_epsilon                     = 0
    # :UserLimit    
    time_limit                      = Inf  
    mip_gap                         = 1e-4
    best_obj_stop                   = NaN
    solution_limit                  = 0
    all_solutions                   = false
    list_of_solutions               = false
    # Parallel  
    processors                      = 1
    # Traversing    
    traverse_strategy               = :BFS
    return SolverOptions(log_levels,branch_strategy,gain_mu,strong_branching_perc,strong_branching_nsteps,strong_restart,
        reliability_branching_threshold, reliability_branching_perc, incumbent_constr,obj_epsilon,time_limit,mip_gap,best_obj_stop,solution_limit,all_solutions,
        list_of_solutions,processors,traverse_strategy)
end

function combine_options(options)
    branch_strategies = Dict{Symbol,Bool}()
    for strat in [:StrongPseudoCost,:PseudoCost,:Reliability,:MostInfeasible]
        branch_strategies[strat] = true
    end

    traverse_strategies = Dict{Symbol,Bool}()
    traverse_strategies[:BFS] = true
    traverse_strategies[:DFS] = true
    traverse_strategies[:DBFS] = true

    options_dict = Dict{Symbol,Any}()
    for kv in options
        if !in(kv[1], fieldnames(SolverOptions))
            warn("Option "*string(kv[1])*" is not available")
        end
        options_dict[kv[1]] = kv[2]
    end
    if haskey(options_dict, :log_levels)
        if length(options_dict[:log_levels]) == 0
            options_dict[:log_levels] = Symbol[]
        end
    end
    defaults = get_default_options()
    for fname in fieldnames(SolverOptions)
        if haskey(options_dict, fname)
            # check branch strategy
            if fname == :branch_strategy 
                if !haskey(branch_strategies, options_dict[fname])
                    error("Branching strategy "*string(options_dict[fname])* " is not supported")
                end
            end

             # check traverse strategy
             if fname == :traverse_strategy 
                if !haskey(traverse_strategies, options_dict[fname])
                    error("Traverse strategy "*string(options_dict[fname])* " is not supported")
                end
            end

            if fieldtype(SolverOptions, fname) != typeof(options_dict[fname])
                options_dict[fname] = convert(fieldtype(SolverOptions,fname), options_dict[fname])
            end
            setfield!(defaults, fname, options_dict[fname])
        end
    end
    return defaults
end

function MINLPBnBSolver(nl_solver::MathProgBase.AbstractMathProgSolver;options...)
    options_obj = combine_options(options)
    return MINLPBnBSolverObj(nl_solver,options_obj)
end
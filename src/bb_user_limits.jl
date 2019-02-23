function isbreak_mip_gap(tree)
    if isdefined(tree, :incumbent) && !tree.options.all_solutions
        incu = tree.incumbent
        b = tree.best_bound
        f = incu.objval
        gap = abs(b-f)/abs(f)
        if gap <= tree.options.mip_gap
            default_opts = get_default_options()
            if tree.options.mip_gap > default_opts.mip_gap
                tree.incumbent = Incumbent(incu.objval, incu.solution, MOI.OBJECTIVE_LIMIT, tree.best_bound)
            else
                tree.incumbent = Incumbent(incu.objval, incu.solution, MOI.LOCALLY_SOLVED, tree.best_bound)
            end
            return true
        end
    end
    return false
end

"""
    isbreak_new_incumbent_limits(tree)

Return true if mip_gap or best_obj_stop is reached
"""
function isbreak_new_incumbent_limits(tree)
    if !isnan(tree.options.best_obj_stop)
        inc_val = tree.incumbent.objval
        bos = tree.options.best_obj_stop
        sense = tree.m.obj_sense
        if (sense == :Min && inc_val <= bos) || (sense == :Max && inc_val >= bos) 
            incu = tree.incumbent
            # TODO: Check if there is a better limit
            tree.incumbent = Incumbent(incu.objval, incu.solution, MOI.OTHER_LIMIT, tree.best_bound)
            return true
        end
    end

    return isbreak_mip_gap(tree)
end

"""
    isbreak_time_limit!(tree)

Check if time limit is reached and  set or update the Incumbent
"""
function isbreak_time_limit!(tree)
    if !isnan(tree.options.time_limit) && time()-tree.start_time >= tree.options.time_limit
        if !isdefined(tree,:incumbent)
            tree.incumbent = Incumbent(NaN, zeros(tree.m.num_var), MOI.TIME_LIMIT, tree.best_bound)
            return true
        else
            tree.incumbent.status = MOI.TIME_LIMIT
            tree.incumbent.best_bound = tree.best_bound
            return true
        end
    end
    return false
end

"""
    isbreak_after_step!(tree)

Check if break...
Break if 
- solution limit is reached
- time limit is reached    
"""
function isbreak_after_step!(tree)
    # maybe break on solution_limit (can be higher if two solutions found in last step)
    if tree.options.solution_limit > 0 && tree.nsolutions >= tree.options.solution_limit
        tree.incumbent.status = MOI.SOLUTION_LIMIT
        tree.incumbent.best_bound = tree.best_bound
        return true
    end
    if isbreak_time_limit!(tree)
        return true
    end

    return false
end
function isbreak_mip_gap(tree)
    if typeof(tree.incumbent) != Void && !tree.options.all_solutions
        b = tree.best_bound
        f = tree.incumbent.objval
        gap_perc = abs(b-f)/abs(f)*100
        if gap_perc <= tree.options.mip_gap
            incu = tree.incumbent
            if tree.options.mip_gap > 1e-2 # bigger than default 
                tree.incumbent = IncumbentSolution(incu.objval,incu.solution,:UserLimit,tree.best_bound)
            else
                tree.incumbent = IncumbentSolution(incu.objval,incu.solution,:Optimal,tree.best_bound)
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
            tree.incumbent = IncumbentSolution(incu.objval,incu.solution,:UserLimit,tree.best_bound)
            return true
        end
    end

    return isbreak_mip_gap(tree)
end

"""
    isbreak_time_limit!(tree)

Check if time limit is reached and  set or update the IncumbentSolution
"""
function isbreak_time_limit!(tree)
    if !isnan(tree.options.time_limit) && time()-tree.start_time >= tree.options.time_limit
        if tree.incumbent == nothing
            tree.incumbent = IncumbentSolution(NaN,zeros(tree.m.num_var),:UserLimit,tree.best_bound)
            return true
        else
            tree.incumbent = IncumbentSolution(tree.incumbent.objval,tree.incumbent.solution,:UserLimit,tree.best_bound)
            return true
        end
    end
    return false
end

"""
    isbreak_after_step!(tree)

Check if break...
Break if 
- incumbent equals best bound
- solution limit is reached
- time limit is reached    
"""
function isbreak_after_step!(tree)
    if !tree.options.all_solutions && tree.incumbent != nothing && tree.incumbent.objval == tree.best_bound
        return true
    end

    # maybe break on solution_limit (can be higher if two solutions found in last step)
    if tree.options.solution_limit > 0 && tree.nsolutions >= tree.options.solution_limit
        incu = tree.incumbent
        tree.incumbent = IncumbentSolution(incu.objval,incu.solution,:UserLimit,tree.best_bound)
        return true
    end
    if isbreak_time_limit!(tree)
        return true
    end

    return false
end
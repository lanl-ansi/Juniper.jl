function isbreak_mip_gap(tree)
    if isdefined(tree, :incumbent) && !tree.options.all_solutions
        incu = tree.incumbent
        b = tree.best_bound
        f = incu.objval
        gap = abs(b - f) / abs(f)
        # if f and g are 0 the gap would be NaN but it's gap=0
        if isapprox(b, f, atol = tree.options.atol)
            gap = 0 # doesn't really equal 0 but should break 
        end
        if gap <= tree.options.mip_gap
            default_opts = get_default_options()
            if tree.options.mip_gap > default_opts.mip_gap
                tree.limit = :MipGap
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
        if (sense == :Min && inc_val <= bos) ||
           (sense == :Max && inc_val >= bos)
            tree.limit = :BestObjStop
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
    if !isnan(tree.options.time_limit) &&
       time() - tree.start_time >= tree.options.time_limit
        tree.limit = :Time
        return true
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
    if tree.options.solution_limit > 0 &&
       tree.nsolutions >= tree.options.solution_limit
        tree.limit = :EnoughSolutions
        return true
    end
    if isbreak_time_limit!(tree)
        return true
    end

    return false
end

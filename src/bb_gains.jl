"""
    update_gains!(tree::BnBTreeObj,node::BnBNode)

Update the objective gains for the branch variable used for node
"""
function update_gains!(tree::BnBTreeObj,parent::BnBNode,l_nd,r_nd)
    
    gain_l = sigma_minus(parent,l_nd,tree.m.solution[parent.var_idx])
    gain_r = sigma_plus(parent,r_nd,tree.m.solution[parent.var_idx])
    idx = tree.var2int_idx[parent.var_idx]

    gain = 0.0
    gain_c = 0
    if !isinf(gain_l) 
        tree.obj_gain_m[idx] = gain_l
        tree.obj_gain_mc[idx] += 1
        gain += gain_l
        gain_c += 1
    end

    if !isinf(gain_r) 
        tree.obj_gain_p[idx] = gain_r
        tree.obj_gain_pc[idx] += 1
        gain += gain_r
        gain_c += 1
    end
    if gain_c == 0
        return 0
    else
        return gain/gain_c
    end
end

"""
    upd_gains_step!(tree,step_obj)

Update the gains using the step_obj if using StrongPseudoCost or PseudoCost
"""
function upd_gains_step!(tree,step_obj)
    branch_strat = tree.options.branch_strategy
    opts = tree.options
    if branch_strat == :StrongPseudoCost && step_obj.counter <= opts.strong_branching_nsteps
        tree.obj_gain_m += step_obj.gains_m
        tree.obj_gain_mc += step_obj.gains_mc
        tree.obj_gain_p += step_obj.gains_p
        tree.obj_gain_pc += step_obj.gains_pc
        strong_int_vars = step_obj.strong_int_vars
        if step_obj.counter == 1
            # all other variables that haven't been checked get the median value of the others
            med_gain_m = median(tree.obj_gain_m[strong_int_vars])
            med_gain_p = median(tree.obj_gain_p[strong_int_vars])
            rest = filter(i->!(i in strong_int_vars),1:tree.m.num_int_bin_var)
            tree.obj_gain_m[rest] += med_gain_m
            tree.obj_gain_p[rest] += med_gain_p
            tree.obj_gain_mc[rest] += 1
            tree.obj_gain_pc[rest] += 1
        end
    elseif branch_strat == :PseudoCost || (branch_strat == :StrongPseudoCost && step_obj.counter > opts.strong_branching_nsteps)
        upd_start = time()
        step_obj.gain_gap = update_gains!(tree,step_obj.node,step_obj.l_nd,step_obj.r_nd)    
        step_obj.upd_gains_time = time()-upd_start
    end
end
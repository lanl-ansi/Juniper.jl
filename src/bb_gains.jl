"""
    compute_gain(node;l_nd::BnBNode=node.left,r_nd::BnBNode=node.right)

Compute the gain of the children to it's parent.
gain = abs(ob-nb)/abs(frac_val-int_val) where ob = old best bound, nb = new best bound
If the state of one child is Integral or Infeasible
return Inf
else return the smaller gain of both children
"""
function compute_gain(node;l_nd::BnBNode=node.left,r_nd::BnBNode=node.right,inf=false)
    gc = 0
    gain_l = 0.0
    gain_r = 0.0
    frac_val = node.solution[node.var_idx]
    if inf && (l_nd.state == :Integral || r_nd.state == :Integral || l_nd.relaxation_state == :Optimal || r_nd.relaxation_state == :Optimal)
        return Inf
    end
    if l_nd.state == :Error && r_nd.state == :Error
        return 0.0
    end
    if l_nd.state == :Branch || l_nd.state == :Integral
        int_val = floor(frac_val)
        gain_l = abs(node.best_bound-l_nd.best_bound)/abs(frac_val-int_val) 
        gc += 1
    end
    if r_nd.state == :Branch || r_nd.state == :Integral
        int_val = ceil(frac_val)
        gain_r = abs(node.best_bound-r_nd.best_bound)/abs(frac_val-int_val)
        gc += 1
    end
    gc == 0 && return 0.0
    # use always minimum of both
    return gain_r > gain_l ? gain_l : gain_r
end

"""
    update_gains!(tree::BnBTreeObj,node::BnBNode,counter)

Update the objective gains for the branch variable used for node
"""
function update_gains!(tree::BnBTreeObj,parent::BnBNode,l_nd,r_nd,counter)
    gain = compute_gain(parent;l_nd=l_nd,r_nd=r_nd)
    
    idx = tree.var2int_idx[parent.var_idx]
    guess = tree.obj_gain[idx]/tree.obj_gain_c[idx]
    if gain == 0 && guess == 0
        gap = 0.0
    elseif gain == 0
        gap = Inf
    else
        gap = abs(guess-gain)/gain*100    
    end

    # update all (just average of the one branch we have)
    if counter == 1
        tree.obj_gain += gain
    else
        tree.obj_gain[idx] += gain
        tree.obj_gain_c[idx] += 1
    end
    return gap
end

"""
    upd_gains_step!(tree,step_obj)

Update the gains using the step_obj if using StrongPseudoCost or PseudoCost
"""
function upd_gains_step!(tree,step_obj)
    branch_strat = tree.options.branch_strategy
    opts = tree.options
    if branch_strat == :StrongPseudoCost && step_obj.counter <= opts.strong_branching_nsteps
        tree.obj_gain += step_obj.gains
        strong_int_vars = step_obj.strong_int_vars
        if step_obj.counter == 1
            # all other variables that haven't been checked get the median value of the others
            med_gain = median(tree.obj_gain[strong_int_vars])
            rest = filter(i->!(i in strong_int_vars),1:tree.m.num_int_bin_var)
            tree.obj_gain[rest] += med_gain
            tree.obj_gain_c += 1
        else
            tree.obj_gain_c[strong_int_vars] += 1
        end
    elseif branch_strat == :PseudoCost || (branch_strat == :StrongPseudoCost && step_obj.counter > opts.strong_branching_nsteps)
        upd_start = time()
        step_obj.gain_gap = update_gains!(tree,step_obj.node,step_obj.l_nd,step_obj.r_nd,step_obj.counter)    
        step_obj.upd_gains_time = time()-upd_start
    end
end
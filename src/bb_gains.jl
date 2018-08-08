importall Base.Operators
function Base.:+(a::GainObj, b::GainObj)
    new_minus = a.minus + b.minus
    new_plus = a.plus + b.plus 
    new_minus_counter = a.minus_counter + b.minus_counter 
    new_plus_counter = a.plus_counter + b.plus_counter 
    return GainObj(new_minus, new_plus, new_minus_counter, new_plus_counter)
end

function init_gains(num_disc_var)
    gains_m = zeros(num_disc_var)
    gains_p = zeros(num_disc_var)
    gains_mc = zeros(Int64,num_disc_var)
    gains_pc = zeros(Int64,num_disc_var)
    return GainObj(gains_m, gains_p, gains_mc, gains_pc)
end

"""
    update_gains!(tree::BnBTreeObj, parent::BnBNode, l_nd, r_nd)

Update the objective gains for the branch variable used for node
"""
function update_gains!(tree::BnBTreeObj, parent::BnBNode, l_nd, r_nd)
    gain_l = sigma_minus(parent, l_nd, parent.solution[parent.var_idx])
    gain_r = sigma_plus(parent,  r_nd, parent.solution[parent.var_idx])
    idx = tree.var2disc_idx[parent.var_idx]

    gain = 0.0
    gain_c = 0
    if !isinf(gain_l) 
        tree.obj_gain.minus[idx] = gain_l
        tree.obj_gain.minus_counter[idx] += 1
        gain += gain_l
        gain_c += 1
    end

    if !isinf(gain_r) 
        tree.obj_gain.plus[idx] = gain_r
        tree.obj_gain.plus_counter[idx] += 1
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
    upd_gains_step!(tree, step_obj)

Update the gains using the step_obj if using StrongPseudoCost or PseudoCost
"""
function upd_gains_step!(tree, step_obj)
    branch_strat = tree.options.branch_strategy
    opts = tree.options
    if step_obj.upd_gains == :GainsToTree || (branch_strat == :StrongPseudoCost && step_obj.counter <= opts.strong_branching_nsteps)
        tree.obj_gain += step_obj.obj_gain
        if step_obj.counter == 1
            cum_counter = tree.obj_gain.minus_counter .+ tree.obj_gain.plus_counter
            strong_int_vars = find(cum_counter .> 0)
            step_obj.strong_int_vars = strong_int_vars
            # all other variables that haven't been checked get the median value of the others
            med_gain_m = median(tree.obj_gain.minus[strong_int_vars])
            med_gain_p = median(tree.obj_gain.plus[strong_int_vars])
            rest = filter(i->!(i in strong_int_vars),1:tree.m.num_disc_var)
            tree.obj_gain.minus[rest] = med_gain_m
            tree.obj_gain.plus[rest] = med_gain_p
            tree.obj_gain.minus_counter[rest] = 1
            tree.obj_gain.plus_counter[rest] = 1
        end
    elseif step_obj.upd_gains == :GuessAndUpdate || branch_strat == :PseudoCost || (branch_strat == :StrongPseudoCost && step_obj.counter > opts.strong_branching_nsteps)
        upd_start = time()
        guess_gain_val = guess_gain(tree, step_obj)
        gain = update_gains!(tree, step_obj.node, step_obj.l_nd, step_obj.r_nd)    
        step_obj.gain_gap = abs(guess_gain_val-gain)/abs(gain)
        if isnan(step_obj.gain_gap) # 0/0
            step_obj.gain_gap = 0.0
        end
        step_obj.upd_gains_time = time()-upd_start
    end
end

function guess_gain(tree, step_obj)
    i = step_obj.var_idx
    inti = tree.var2disc_idx[i]
    x = step_obj.node.solution
    g_minus, g_minus_c = tree.obj_gain.minus,tree.obj_gain.minus_counter
    g_plus, g_plus_c = tree.obj_gain.plus,tree.obj_gain.plus_counter
    g_minus_c += map(i -> (i == 0) && (i = 1), g_minus_c)
    g_plus_c += map(i -> (i == 0) && (i = 1), g_plus_c)
    mu = tree.options.gain_mu
    return score(f_minus(x[i])*g_minus[inti]/g_minus_c[inti],f_plus(x[i])*g_plus[inti]/g_plus_c[inti],mu)
end
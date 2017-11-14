"""
branch_mostinfeasible(num_var,var_type,node)

Get the index of an integer variable which is currently continuous which is most unintegral.
(nearest to *.5)
"""
function branch_mostinfeasible(m,node)
    x = node.solution
    idx = 0
    max_diff = 0
    for i=1:m.num_var
        if m.var_type[i] != :Cont
            diff = abs(x[i]-round(x[i]))
            if diff > max_diff
                idx = i
                max_diff = diff
            end
        end
    end
    return idx
end

"""
init_strong_restart!(node, var_idx, int_var_idx, l_nd, r_nd, reasonable_int_vars, infeasible_int_vars)

Tighten the bounds for the node and check if there are variables that need to be checked for a restart.
"""
function init_strong_restart!(node, var_idx, int_var_idx, l_nd, r_nd, reasonable_int_vars, infeasible_int_vars)
    restart = false

    # set the bounds directly for the node
    # also update the best bound and the solution
    if l_nd.relaxation_state != :Optimal
        node.l_var[var_idx] = ceil(node.solution[var_idx])
        node.best_bound = r_nd.best_bound
        node.solution = r_nd.solution
    else
        node.u_var[var_idx] = floor(node.solution[var_idx])
        node.best_bound = l_nd.best_bound
        node.solution = l_nd.solution
    end

    push!(infeasible_int_vars,int_var_idx)

    if length(reasonable_int_vars) == length(infeasible_int_vars)
        # basically branching on the last infeasible variable 
        max_gain_var = infeasible_int_vars[end]
        strong_int_vars = [infeasible_int_vars[end]] # don't divide by 0 later
    else
        max_gain_var = 0
        strong_int_vars = zeros(Int64,0)
        restart = true
    end
    return restart, infeasible_int_vars, max_gain_var, strong_int_vars
end

function branch_strong_on(m,opts,step_obj,
        reasonable_int_vars, int2var_idx, strong_restart, counter)

    max_gain = -Inf # then one is definitely better
    max_gain_var = 0
    strong_int_vars = zeros(Int64,0)        
    strong_restarts = -1
    restart = true
    status = :Normal
    node = step_obj.node
    infeasible_int_vars = zeros(Int64,0)
    gains_m = zeros(m.num_int_bin_var)
    gains_p = zeros(m.num_int_bin_var)
    gains_mc = zeros(Int64,m.num_int_bin_var)
    gains_pc = zeros(Int64,m.num_int_bin_var)

    left_node, right_node = nothing, nothing
    while restart 
        strong_restarts += 1 # is init with -1
        restart = false
        for int_var_idx in reasonable_int_vars
            # don't rerun if the variable has already one infeasible node
            if int_var_idx in infeasible_int_vars
                continue
            end
            push!(strong_int_vars, int_var_idx)
            var_idx = int2var_idx[int_var_idx]
            step_obj.var_idx = var_idx
            u_b, l_b = node.u_var[var_idx], node.l_var[var_idx]
            # don't rerun if bounds are exact or is type correct
            if isapprox(u_b,l_b,atol=atol) || is_type_correct(node.solution[var_idx],m.var_type[var_idx])
                continue
            end
            # branch on the current variable and get the corresponding children
            l_nd,r_nd = branch!(m,opts,step_obj,counter;temp=true)
            if l_nd.relaxation_state != :Optimal && r_nd.relaxation_state != :Optimal && counter == 1
                # TODO: Might be Error instead of infeasible
                status = :GlobalInfeasible
                break
            end

            # if restart is true => check if one part is infeasible => update bounds & restart
            if strong_restart == true
                if l_nd.relaxation_state != :Optimal || r_nd.relaxation_state != :Optimal
                    max_gain = 0.0
                    if l_nd.relaxation_state != :Optimal && r_nd.relaxation_state != :Optimal
                        # TODO: Might be Error instead of infeasible
                        status = :LocalInfeasible
                        break
                    end
                    restart,infeasible_int_vars,max_gain_var,strong_int_vars = init_strong_restart!(node, var_idx, int_var_idx, l_nd, r_nd, reasonable_int_vars, infeasible_int_vars)
                    break
                end
            end
            println("x: ", node.solution[var_idx])
            gain_l = sigma_minus(node,l_nd,m.solution[node.var_idx])
            println("gain_l: ", gain_l)
            gain_r = sigma_plus(node,r_nd,m.solution[node.var_idx])
            println("gain_r: ", gain_r)
            gain = (gain_l+gain_r)/2
            if gain > max_gain
                max_gain = gain
                max_gain_var = var_idx
                left_node = l_nd
                right_node = r_nd
                # gain is set to inf if Integral or Infeasible
                # TODO: Might be reasonable to use something different
                if gain == Inf
                    break
                end
            end
            if !isinf(gain_l)
                gains_m[int_var_idx] = gain_l
                gains_mc[int_var_idx] = 1
            end

            if !isinf(gain_r)
                gains_p[int_var_idx] = gain_r
                gains_pc[int_var_idx] = 1
            end
        end
    end
    return status, max_gain_var, left_node, right_node, 
    (gains_m, gains_mc, gains_p, gains_pc), strong_restarts, strong_int_vars
end

"""
branch_strong!(m,opts,int2var_idx,step_obj,counter)

Try to branch on a few different variables and choose the one with highest obj_gain.
Update obj_gain for the variables tried and average the other ones.
"""
function branch_strong!(m,opts,int2var_idx,step_obj,counter)
    node = step_obj.node

    # generate an of variables to branch on
    num_strong_var = Int(round((opts.strong_branching_perc/100)*m.num_int_bin_var))
    # if smaller than 2 it doesn't make sense
    num_strong_var = num_strong_var < 2 ? 2 : num_strong_var

    # get reasonable candidates (not type correct and not already perfectly bounded)
    int_vars = m.num_int_bin_var
    reasonable_int_vars = zeros(Int64,0)
    for i=1:int_vars
        idx = int2var_idx[i]
        u_b = node.u_var[idx]
        l_b = node.l_var[idx]
        if isapprox(u_b,l_b,atol=atol) || is_type_correct(node.solution[idx],m.var_type[idx])
            continue
        end
        push!(reasonable_int_vars,i)
    end
    shuffle!(reasonable_int_vars)
    reasonable_int_vars = reasonable_int_vars[1:minimum([num_strong_var,length(reasonable_int_vars)])]

    # compute the gain for each reasonable candidate and choose the highest
    left_node = nothing
    right_node = nothing

    status, max_gain_var,  left_node, right_node, gains, strong_restarts, strong_int_vars = branch_strong_on(m,opts,step_obj,
        reasonable_int_vars, int2var_idx, opts.strong_restart, counter)

    gains_m, gains_mc, gains_p, gains_pc = gains
    step_obj.gains_m += gains_m
    step_obj.gains_mc += gains_mc
    step_obj.gains_p += gains_p
    step_obj.gains_pc += gains_pc

    if status != :GlobalInfeasible && status != :LocalInfeasible
        step_obj.l_nd = left_node
        step_obj.r_nd = right_node
    
        # set the variable to branch (best gain)
        node.state = :Done
        node.var_idx = max_gain_var
        step_obj.strong_int_vars = strong_int_vars
    end

    @assert max_gain_var != 0 || status == :LocalInfeasible || status == :GlobalInfeasible || node.state == :Infeasible
    return status, max_gain_var, strong_restarts
end

function branch_reliable!(m,opts,step_obj,int2var_idx,g_minus,g_minus_c,g_plus,g_plus_c,mu,counter) 
    idx = 0
    node = step_obj.node
    reliability_param = 5
    reliability_perc = 50
    num_strong_var = Int(round((reliability_perc/100)*m.num_int_bin_var))
    # if smaller than 2 it doesn't make sense
    num_strong_var = num_strong_var < 2 ? 2 : num_strong_var

    gmc_r = g_minus_c .< reliability_param
    gpc_r = g_plus_c  .< reliability_param
    
    reasonable_int_vars = []
    for i=1:length(gmc_r)
        if gmc_r[i] || gpc_r[i]
            push!(reasonable_int_vars,i)
        end
    end
    if length(reasonable_int_vars) > 0
        shuffle!(reasonable_int_vars)
        num_reasonable = num_strong_var < length(reasonable_int_vars) ? num_strong_var : length(reasonable_int_vars)
        reasonable_int_vars = reasonable_int_vars[1:num_reasonable]
        
        println("reasonable_int_vars: ", reasonable_int_vars)
        @assert are_type_correct(node.solution,m.var_type) == false
        status, max_gain_var,  left_node, right_node, gains, strong_restarts, strong_int_vars = branch_strong_on(m,opts,step_obj,
            reasonable_int_vars, int2var_idx, opts.strong_restart, counter)
        gains_m, gains_mc, gains_p, gains_pc = gains
        println("gains: ", gains)
        step_obj.gains_m += gains_m
        step_obj.gains_mc += gains_mc
        step_obj.gains_p += gains_p
        step_obj.gains_pc += gains_pc
        step_obj.strong_int_vars = strong_int_vars
        step_obj.upd_gains = :GainsToTree
        new_gains_m = step_obj.gains_m
        new_gains_mc = step_obj.gains_mc
        new_gains_p = step_obj.gains_p
        new_gains_pc = step_obj.gains_pc
    else 
        step_obj.upd_gains = :GuessAndUpdate
        new_gains_m = g_minus
        new_gains_mc = g_minus_c
        new_gains_p = g_plus
        new_gains_pc = g_plus_c
    end
    println("step_obj.upd_gains: ", step_obj.upd_gains)
    scores, sort_idx = sorted_score_idx(m.solution,new_gains_m,new_gains_mc,new_gains_p,new_gains_pc,int2var_idx,mu)
    println("sorted scores: ", scores[sort_idx])
    println("new_gains_m: ", new_gains_m)
    println("new_gains_mc: ", new_gains_mc)
    println("new_gains_p: ", new_gains_p)
    println("new_gains_pc: ", new_gains_pc)
    @assert are_type_correct(node.solution,m.var_type) == false
    for l_idx in sort_idx
        var_idx = int2var_idx[l_idx]
        if !is_type_correct(node.solution[var_idx],m.var_type[var_idx])
            u_b = node.u_var[var_idx]
            l_b = node.l_var[var_idx]
            # if the upper bound is the lower bound => no reason to branch
            if isapprox(u_b,l_b,atol=atol)
                continue
            end
            idx = var_idx
            break
        end
    end
    @assert idx != 0
    return idx
end

function branch_pseudo(m,node,int2var_idx,g_minus,g_minus_c,g_plus,g_plus_c,mu)
    # use the one with highest obj_gain which is currently continous
    idx = 0
    sort_idx = sorted_score_idx(m.solution,g_minus,g_minus_c,g_plus,g_plus_c,int2var_idx,mu)
    for l_idx in sort_idx
        var_idx = int2var_idx[l_idx]
        if !is_type_correct(node.solution[var_idx],m.var_type[var_idx])
            u_b = node.u_var[var_idx]
            l_b = node.l_var[var_idx]
            # if the upper bound is the lower bound => no reason to branch
            if isapprox(u_b,l_b,atol=atol)
                continue
            end
            idx = var_idx
            break
        end
    end
    return idx
end

function sorted_score_idx(x,g_minus,g_minus_c,g_plus,g_plus_c,i2v,mu)
    scores = [score(f_minus(x[i2v[i]])*g_minus[i]/g_minus_c[i],f_plus(x[i2v[i]])*g_plus[i]/g_plus_c[i],mu) for i=1:length(g_minus)]
    sortedidx = sortperm(scores; rev=true)
    return scores,sortedidx
end

"""
    Score function from 
    http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.92.7117&rep=rep1&type=pdf
"""
function score(q_m,q_p,mu) 
    minq = q_m < q_p ? q_m : q_p
    maxq = q_m < q_p ? q_p : q_m
    return (1-mu)*minq+mu*maxq
end

function diff_obj(node,leaf)
    if leaf.relaxation_state == :Optimal
        return abs(node.best_bound - leaf.best_bound)
    else
        return Inf
    end
end

f_plus(x) = ceil(x)-x
f_minus(x) = x-floor(x)
sigma_plus(node,r_nd,x) = diff_obj(node,r_nd)/f_plus(x)
sigma_minus(node,l_nd,x) = diff_obj(node,l_nd)/f_minus(x)
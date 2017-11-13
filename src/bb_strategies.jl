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

"""
branch_strong!(m,opts,int2var_idx,step_obj,counter)

Try to branch on a few different variables and choose the one with highest obj_gain.
Update obj_gain for the variables tried and average the other ones.
"""
function branch_strong!(m,opts,int2var_idx,step_obj,counter)
    function init_variables()
        max_gain = -Inf # then one is definitely better
        max_gain_var = 0
        strong_int_vars = zeros(Int64,0)
        return max_gain, max_gain_var, strong_int_vars
    end

    node = step_obj.node

    strong_restarts = -1 

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
    max_gain, max_gain_var, strong_int_vars = init_variables()
    left_node = nothing
    right_node = nothing

    restart = true
    infeasible_int_vars = zeros(Int64,0)
    status = :Normal
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
            if opts.strong_restart == true
                if l_nd.relaxation_state != :Optimal || r_nd.relaxation_state != :Optimal
                    max_gain = 0.0
                    if l_nd.relaxation_state != :Optimal && r_nd.relaxation_state != :Optimal
                        # TODO: Might be Error instead of infeasible
                        status = :LocalInfeasible
                        break
                    end
                    restart,infeasible_int_vars,max_gain_var,strong_int_vars = init_strong_restart!(node, var_idx, int_var_idx, l_nd, r_nd, reasonable_int_vars, infeasible_int_vars)
                end
            end
            gain_l = sigma_minus(node,l_nd,m.solution[node.var_idx])
            gain_r = sigma_plus(node,r_nd,m.solution[node.var_idx])
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
                step_obj.gains_m[int_var_idx] = gain_l
                step_obj.gains_mc[int_var_idx] += 1
            end

            if !isinf(gain_r)
                step_obj.gains_p[int_var_idx] = gain_r
                step_obj.gains_pc[int_var_idx] += 1
            end
        end
    end

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

function branch_pseudo(m,node,int2var_idx,g_minus,g_minus_c,g_plus,g_plus_c)
    # use the one with highest obj_gain which is currently continous
    idx = 0
    sort_idx = sorted_score_idx(m.solution,g_minus,g_minus_c,g_plus,g_plus_c,int2var_idx)
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

function sorted_score_idx(x,g_minus,g_minus_c,g_plus,g_plus_c,i2v)
    scores = [score(f_minus(x[i2v[i]])*g_minus[i]/g_minus_c[i],f_plus(x[i2v[i]])*g_plus[i]/g_plus_c[i]) for i=1:length(g_minus)]
    sortedidx = sortperm(scores; rev=true)
    return sortedidx
end

"""
    Score function from 
    http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.92.7117&rep=rep1&type=pdf
"""
function score(q_m,q_p) 
    mu = 0.167 # TODO changeable
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
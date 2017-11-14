function init(start_time, m)
    srand(1)
    node = BnBNode(1,1,m.l_var,m.u_var,m.solution,0,:Branch,:Optimal,m.objval)
    obj_gain_m = zeros(m.num_int_bin_var)
    obj_gain_p = zeros(m.num_int_bin_var)
    obj_gain_mc = zeros(Int64,m.num_int_bin_var)
    obj_gain_pc = zeros(Int64,m.num_int_bin_var)
    int2var_idx = zeros(m.num_int_bin_var)
    var2int_idx = zeros(m.num_var)
    int_i = 1
    for i=1:m.num_var
        if m.var_type[i] != :Cont
            int2var_idx[int_i] = i
            var2int_idx[i] = int_i
            int_i += 1
        end
    end
    factor = 1
    if m.obj_sense == :Min
        factor = -1
    end
    return BnBTreeObj(m,nothing,obj_gain_m,obj_gain_p,obj_gain_mc,obj_gain_pc,int2var_idx,var2int_idx,m.options,
                    factor,start_time,0,[node],NaN)
end

function new_default_node(idx,level,l_var,u_var,solution;
                            var_idx=0,
                            state=:Solve,relaxation_state=:Solve,best_bound=nothing)

    l_var = copy(l_var)
    u_var = copy(u_var)
    solution = copy(solution)
    return BnBNode(idx,level,l_var,u_var,solution,var_idx,state,relaxation_state,best_bound)     
end

function new_default_step_obj(m,node)
    gains_m = zeros(m.num_int_bin_var)
    gains_mc = ones(Int64,m.num_int_bin_var)
    gains_p = zeros(m.num_int_bin_var)
    gains_pc = ones(Int64,m.num_int_bin_var)
    idx_time = 0.0
    leaf_idx_time = 0.0
    upd_gains_time = 0.0
    leaf_branch_time = 0.0
    branch_time = 0.0
    return StepObj(node,0,:None,0,0.0,gains_m,gains_mc,gains_p,gains_pc,zeros(Int64,0),
    idx_time,leaf_idx_time,upd_gains_time,leaf_branch_time,branch_time,[],[],nothing,nothing,0,:None)
end
function init(start_time, m)
    srand(1)
    node = BnBNode(1,1,m.l_var,m.u_var,m.solution,0,:Branch,m.objval)
    obj_gain = zeros(m.num_int_bin_var)
    obj_gain_c = zeros(m.num_int_bin_var)
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
    return BnBTreeObj(m,nothing,obj_gain,obj_gain_c,int2var_idx,var2int_idx,m.options,
                    factor,start_time,0,[node],NaN)
end

function new_default_node(idx,level,l_var,u_var,solution;
                            var_idx=0,
                            state=:Solve,best_bound=nothing)

    l_var = copy(l_var)
    u_var = copy(u_var)
    solution = copy(solution)
    return BnBNode(idx,level,l_var,u_var,solution,var_idx,state,best_bound)     
end

function new_default_step_obj(m,node)
    idx_time = 0.0
    leaf_idx_time = 0.0
    upd_gains_time = 0.0
    leaf_branch_time = 0.0
    branch_time = 0.0
    return StepObj(node,0,:None,0,0.0,zeros(m.num_int_bin_var),zeros(Int64,0),idx_time,leaf_idx_time,upd_gains_time,leaf_branch_time,branch_time,[],[],nothing,nothing,0)
end
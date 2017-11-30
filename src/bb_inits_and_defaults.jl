function init(start_time, m; inc_sol = nothing, inc_obj = nothing)
    srand(1)
    node = BnBNode(1, 1, m.l_var, m.u_var, m.solution, 0, :Branch, :Optimal, m.objval)
    obj_gain_m = zeros(m.num_int_bin_var)
    obj_gain_p = zeros(m.num_int_bin_var)
    obj_gain_mc = zeros(Int64, m.num_int_bin_var)
    obj_gain_pc = zeros(Int64, m.num_int_bin_var)
    obj_gain = GainObj(obj_gain_m, obj_gain_p, obj_gain_mc, obj_gain_pc)
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
    bnbTree = BnBTreeObj()
    bnbTree.m           = m
    bnbTree.obj_gain    = obj_gain
    bnbTree.int2var_idx = int2var_idx
    bnbTree.var2int_idx = var2int_idx
    bnbTree.options     = m.options
    bnbTree.obj_fac     = factor
    bnbTree.start_time  = start_time
    bnbTree.nsolutions  = 0
    bnbTree.branch_nodes = [node]
    bnbTree.best_bound  = NaN
    bnbTree.mutex_get_node = false

    if inc_sol != nothing
        bnbTree.incumbent = Incumbent(inc_obj, inc_sol, :Optimal, m.objval)
        
        if m.options.incumbent_constr
            add_incumbent_constr(m,bnbTree.incumbent)
            bnbTree.nsolutions += 1
            m.ncuts += 1
        end
    end

    return bnbTree
end

function new_default_node(idx, level, l_var, u_var, solution;
                            var_idx=0, state=:Solve, relaxation_state=:Solve, best_bound=NaN)

    l_var = copy(l_var)
    u_var = copy(u_var)
    solution = copy(solution)
    return BnBNode(idx, level, l_var, u_var, solution, var_idx, state, relaxation_state, best_bound)
end

function init_time_obj()
    return TimeObj(0.0,0.0,0.0,0.0,0.0)
end

function new_default_step_obj(m,node)
    gains_m = zeros(m.num_int_bin_var)
    gains_mc = zeros(Int64, m.num_int_bin_var)
    gains_p = zeros(m.num_int_bin_var)
    gains_pc = zeros(Int64, m.num_int_bin_var)
    gains = GainObj(gains_m, gains_p, gains_mc, gains_pc)
    idx_time = 0.0
    node_idx_time = 0.0
    upd_gains_time = 0.0
    node_branch_time = 0.0
    branch_time = 0.0
    step_obj = StepObj()
    step_obj.node             = node    
    step_obj.var_idx          = 0    
    step_obj.state            = :None   
    step_obj.nrestarts        = 0   
    step_obj.gain_gap         = 0.0   
    step_obj.obj_gain         = gains
    step_obj.strong_int_vars  = zeros(Int64,0)    
    step_obj.idx_time         = idx_time   
    step_obj.node_idx_time    = node_idx_time   
    step_obj.upd_gains_time   = upd_gains_time   
    step_obj.node_branch_time = node_branch_time   
    step_obj.branch_time      = branch_time   
    step_obj.integral         = []   
    step_obj.branch           = []   
    step_obj.counter          = 0   
    step_obj.upd_gains        = :None
    return step_obj
end
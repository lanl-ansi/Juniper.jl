function init(start_time, m; inc_sol = nothing, inc_obj = nothing)
    srand(1)
    hash_val = string(hash(hcat(m.l_var,m.u_var)))
    node = BnBNode(1, 1, m.l_var, m.u_var, m.solution, 0, :Branch, :Optimal, m.objval,[],hash_val)
    obj_gain_m = zeros(m.num_disc_var)
    obj_gain_p = zeros(m.num_disc_var)
    obj_gain_mc = zeros(Int64, m.num_disc_var)
    obj_gain_pc = zeros(Int64, m.num_disc_var)
    obj_gain = GainObj(obj_gain_m, obj_gain_p, obj_gain_mc, obj_gain_pc)
    disc2var_idx = zeros(m.num_disc_var)
    var2disc_idx = zeros(m.num_var)
    int_i = 1
    for i=1:m.num_var
        if m.var_type[i] != :Cont
            disc2var_idx[int_i] = i
            var2disc_idx[i] = int_i
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
    bnbTree.disc2var_idx = disc2var_idx
    bnbTree.var2disc_idx = var2disc_idx
    bnbTree.options     = m.options
    bnbTree.obj_fac     = factor
    bnbTree.start_time  = start_time
    bnbTree.nsolutions  = 0
    bnbTree.branch_nodes = [node]
    bnbTree.best_bound  = m.objval

    if inc_sol != nothing
        bnbTree.incumbent = Incumbent(inc_obj, inc_sol, :Optimal, m.objval)
        bnbTree.nsolutions += 1
        if m.options.incumbent_constr
            add_incumbent_constr(m,bnbTree.incumbent)
            m.ncuts += 1
        end
    end

    return bnbTree
end

function new_default_node(idx, level, l_var, u_var, solution;
                            var_idx=0, state=:Solve, relaxation_state=:Solve, best_bound=NaN, path=[])

    l_var = copy(l_var)
    u_var = copy(u_var)
    solution = copy(solution)
    hash_val = string(hash(hcat(l_var,u_var)))
    return BnBNode(idx, level, l_var, u_var, solution, var_idx, state, relaxation_state, best_bound, path, hash_val)
end

function new_left_node(node, u_var; 
                       var_idx=0, state=:Solve, relaxation_state=:Solve, best_bound=NaN, path=[])
    l_var = copy(node.l_var)
    u_var = copy(u_var)
    solution = NaN*ones(length(node.solution))
    hash_val = string(hash(hcat(l_var,u_var)))
    return BnBNode(node.idx*2, node.level+1, l_var, u_var, solution, var_idx, state, relaxation_state, best_bound, path, hash_val)
end

function new_right_node(node, l_var; 
                       var_idx=0, state=:Solve, relaxation_state=:Solve, best_bound=NaN, path=[])
    l_var = copy(l_var)
    u_var = copy(node.u_var)
    solution = NaN*ones(length(node.solution))
    hash_val = string(hash(hcat(l_var,u_var)))
    return BnBNode(node.idx*2+1, node.level+1, l_var, u_var, solution, var_idx, state, relaxation_state, best_bound, path, hash_val)
end

function init_time_obj()
    return TimeObj(0.0,0.0,0.0,0.0,0.0)
end

function new_default_step_obj(m,node)
    gains_m = zeros(m.num_disc_var)
    gains_mc = zeros(Int64, m.num_disc_var)
    gains_p = zeros(m.num_disc_var)
    gains_pc = zeros(Int64, m.num_disc_var)
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
    step_obj.idx_time         = idx_time   
    step_obj.node_idx_time    = node_idx_time   
    step_obj.upd_gains_time   = upd_gains_time   
    step_obj.node_branch_time = node_branch_time   
    step_obj.branch_time      = branch_time   
    step_obj.integral         = []   
    step_obj.branch           = []   
    step_obj.counter          = 0   
    step_obj.upd_gains        = :None
    step_obj.strong_int_vars  = Int64[]
    return step_obj
end
step_obj_primitives = [:var_idx,:state,:nrestarts,:gain_gap,
                       :strong_disc_vars,:idx_time,:node_idx_time,:upd_gains_time,:branch_time,
                       :counter,:upd_gains]
node_primitives = [:level,:var_idx,:l_var,:u_var,:solution,:state,:relaxation_state,:best_bound]
gain_obj_primitives = [:minus,:plus,:minus_counter,:plus_counter]

typedict(x,keys) = Dict(fn=>getfield(x, fn) for fn ∈ keys) 

function get_entry_dict(step_obj)
    node = step_obj.node
    d = Dict{Symbol,Any}()
    d_l = Dict{Symbol,Any}()
    d_r = Dict{Symbol,Any}()
    d[:hash] = node.hash
    d[:children] = Vector{Dict}()
    n = Dict{Symbol,Any}()
    n_l = Dict{Symbol,Any}()
    n_r = Dict{Symbol,Any}()

    step_obj_dict = typedict(step_obj, step_obj_primitives)
    node_dict = typedict(node, node_primitives)
    gain_obj_dict = typedict(step_obj.obj_gain, gain_obj_primitives)

    step_obj_dict[:node] = node_dict
    step_obj_dict[:obj_gain] = gain_obj_dict
    d[:step_obj] = step_obj_dict

    d_l[:hash] = step_obj.l_nd.hash
    n_l[:state] = step_obj.l_nd.state
    n_l[:best_bound] = step_obj.l_nd.best_bound
    n_l[:relaxation_state] = step_obj.l_nd.relaxation_state
    if n_l[:state] == :Integral
        n_l[:solution] = step_obj.l_nd.solution
    end
    d_r[:hash] = step_obj.r_nd.hash
    n_r[:state] = step_obj.r_nd.state
    n_r[:best_bound] = step_obj.r_nd.best_bound
    n_r[:relaxation_state] = step_obj.r_nd.relaxation_state
    if n_r[:state] == :Integral
        n_r[:solution] = step_obj.r_nd.solution
    end
    d_l[:step_obj] = Dict{Symbol,Any}()
    d_r[:step_obj] = Dict{Symbol,Any}()
    d_l[:step_obj][:node] = n_l
    d_r[:step_obj][:node] = n_r
    push!(d[:children],d_l)
    push!(d[:children],d_r)
    return d
end

function upd_node_dict!(cd, step_obj)
    ed = get_entry_dict(step_obj)
    cd[:children] = ed[:children]
    cd[:step_obj] = ed[:step_obj]
    cd[:hash] = ed[:hash]
end


function push_step2treeDict!(d, step_obj)
    c = step_obj.counter
    node = step_obj.node
    if length(node.path) == 0
        d = get_entry_dict(step_obj)
    else 
        path = copy(node.path)
        cd = d
        pnode = popfirst!(path)
        push!(path,step_obj.node)
        while length(path) > 0
            pnode = popfirst!(path)
            phash = pnode.hash
            if cd[:children][1][:hash] == phash
                cd = cd[:children][1]
            else 
                cd = cd[:children][2]
            end
        end
        upd_node_dict!(cd,step_obj)
    end
    return d
end

function debug_init(d)
    d[:relaxation] = Dict{Symbol,Any}()
    d[:info] = Dict{Symbol,Any}()
end 

function debug_fill_basic(d,m,restarts)
    d[:relaxation][:status] = m.status
    d[:relaxation][:time] = m.relaxation_time
    d[:relaxation][:nrestarts] = restarts
    d[:info][:sense] = m.obj_sense
    d[:info][:nintvars] = m.nintvars
    d[:info][:nbinvars] = m.nbinvars
    d[:info][:var2disc_idx] = m.var2disc_idx
end


function debug_objective(d,m)
    d[:relaxation][:objval] = m.objval
    d[:relaxation][:solution] = m.solution
end 
    
function debug_restart_values(d,restart_vals)
    if !haskey(d[:relaxation], :restarts)
        d[:relaxation][:restarts] = []
    end
    push!(d[:relaxation][:restarts], restart_vals)
end

function debug_set_solution(d,m)
    d[:solution] = Dict{Symbol,Any}()
    d[:solution][:objval] = m.objval
    d[:solution][:best_bound] = m.best_bound
    d[:solution][:status] = m.status
    d[:solution][:solution] = m.solution
    d[:solution][:time] = m.soltime
end

function debug_set_tree_obj_gain!(tree::BnBTreeObj)
    tree.m.debugDict[:obj_gain] = zeros(4,tree.m.num_disc_var)
    tree.m.debugDict[:obj_gain][1,:] = tree.obj_gain.minus
    tree.m.debugDict[:obj_gain][2,:] = tree.obj_gain.plus
    tree.m.debugDict[:obj_gain][3,:] = tree.obj_gain.minus_counter
    tree.m.debugDict[:obj_gain][4,:] = tree.obj_gain.plus_counter
end
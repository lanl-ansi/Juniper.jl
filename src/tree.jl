function get_node_dict(step_obj)
    node = step_obj.node
    d = Dict{Symbol,Any}()
    d_l = Dict{Symbol,Any}()
    d_r = Dict{Symbol,Any}()
    d[:hash] = node.hash
    d[:children] = Vector{Dict}()
    n = Dict{Symbol,Any}()
    n_l = Dict{Symbol,Any}()
    n_r = Dict{Symbol,Any}()

    n[:counter] = step_obj.counter
    n[:best_bound] = node.best_bound
    n[:var_idx] = node.var_idx
    n[:state] = node.state
    d[:node] = n
    d_l[:hash] = step_obj.l_nd.hash
    n_l[:state] =step_obj.l_nd.state
    n_l[:best_bound] =step_obj.l_nd.best_bound
    n_l[:rel_state] =step_obj.l_nd.relaxation_state
    d_r[:hash] = step_obj.r_nd.hash
    n_r[:state] =step_obj.r_nd.state
    n_r[:best_bound] =step_obj.r_nd.best_bound
    n_r[:rel_state] =step_obj.r_nd.relaxation_state
    d_l[:node] = n_l
    d_r[:node] = n_r
    push!(d[:children],d_l)
    push!(d[:children],d_r)
    return d
end

function upd_node_dict!(cd, step_obj)
    node_dict = get_node_dict(step_obj)
    cd[:node][:var_idx] = node_dict[:node][:var_idx]
    cd[:node][:state] = node_dict[:node][:state]
    cd[:children] = node_dict[:children]
    cd[:node][:counter] = node_dict[:node][:counter]
end


function push_step2treeDict!(d, step_obj)
    c = step_obj.counter
    node = step_obj.node
    if length(node.path) == 0
        d = get_node_dict(step_obj)
    else 
        path = copy(node.path)
        cd = d
        pnode = shift!(path)
        push!(path,step_obj.node)
        while length(path) > 0
            pnode = shift!(path)
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

function debug_init(d,m,restarts)
    d[:relaxation] = Dict{Symbol,Any}()
    d[:relaxation][:status] = m.status
    d[:relaxation][:time] = m.relaxation_time
end 

function debug_objective(d,m)
    d[:relaxation][:objval] = m.objval
    d[:relaxation][:solution] = m.solution
end 
    
function get_node_dict(step_obj)
    node = step_obj.node
    d = Dict{String,Any}()
    d["hash"] = node.hash
    d["counter"] = step_obj.counter
    d["best_bound"] = node.best_bound
    d["var_idx"] = node.var_idx
    d["state"] = node.state
    d["children"] = Vector{Dict}()
    d_l = Dict{Any,Any}()
    d_l["hash"] = step_obj.l_nd.hash
    d_l["state"] =step_obj.l_nd.state
    d_l["best_bound"] =step_obj.l_nd.best_bound
    d_l["rel_state"] =step_obj.l_nd.relaxation_state
    d_r = Dict{Any,Any}()
    d_r["hash"] = step_obj.r_nd.hash
    d_r["state"] =step_obj.r_nd.state
    d_r["best_bound"] =step_obj.r_nd.best_bound
    d_r["rel_state"] =step_obj.r_nd.relaxation_state
    push!(d["children"],d_l)
    push!(d["children"],d_r)
    return d
end

function upd_node_dict!(cd, step_obj)
    node_dict = get_node_dict(step_obj)
    cd["var_idx"] = node_dict["var_idx"]
    cd["state"] = node_dict["state"]
    cd["children"] = node_dict["children"]
    cd["counter"] = node_dict["counter"]
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
            if cd["children"][1]["hash"] == phash
                cd = cd["children"][1]
            else 
                cd = cd["children"][2]
            end
        end
        upd_node_dict!(cd,step_obj)
    end
    return d
end
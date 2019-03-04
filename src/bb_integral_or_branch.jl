"""
    push_integral_or_branch!(step_obj, cnode)

Add integral or branch node to step_obj
"""
function push_integral_or_branch!(step_obj, cnode)
    # check if all int vars are int
    if cnode.state == :Integral
        push!(step_obj.integral, cnode)
    else
        push!(step_obj.branch, cnode)
    end
end

"""
    set_cnode_state!(cnode, m, step_obj, disc2var_idx)

Set the state of the current node to :Integral or :Branch
"""
function set_cnode_state!(cnode, m, step_obj, disc2var_idx)
    # check if all int vars are int
    if are_type_correct(cnode.solution, m.var_type, disc2var_idx, m.options.atol)
        if cnode.relaxation_state == MOI.ALMOST_LOCALLY_SOLVED && m.options.allow_almost_solved_integral
            @warn "Integral leaf node only almost locally solved. Disallowable with `allow_almost_solved_integral=false`"
        end
        cnode.state = :Integral
    else
        cnode.state = :Branch
    end
end

"""
    new_integral!(tree, node)

Update the incumbent and add obj constr if in options
Add to solutions if list_of_solutions in options
"""
function new_integral!(tree, node)
    # node.best_bound is the objective for integral values
    tree.nsolutions += 1
    if tree.options.list_of_solutions
        push!(tree.m.solutions, Juniper.SolutionObj(node.solution, node.best_bound))
    end
    if update_incumbent!(tree,node) # returns if new 
        if tree.options.incumbent_constr
            if tree.options.processors > 1
                np = nprocs()  # determine the number of processes available
                if tree.options.processors+1 < np
                    np = tree.options.processors+1
                end
                for p=2:np
                    sendto(p, is_newincumbent=true)
                end
            end
            add_incumbent_constr(tree.m, tree.incumbent)
            tree.m.ncuts += 1
        end
    end
end

"""
    push_to_branch_list!(tree, node)

Push a node to the list of branch nodes if better than incumbent
or if all solutions are wanted
"""
function push_to_branch_list!(tree, node)
    incu = isdefined(tree,:incumbent) ? tree.incumbent : false
    if tree.options.all_solutions || !isdefined(tree,:incumbent) || tree.obj_fac*node.best_bound >= tree.obj_fac*incu.objval
        push!(tree.branch_nodes, node)
    end
end

"""
    upd_integral_branch!(tree, step_obj)

Update the list of integral and branch nodes using the new step_obj
Return true if break
"""
function upd_integral_branch!(tree, step_obj)
    for integral_node in step_obj.integral
        new_integral!(tree,integral_node)
        if isbreak_new_incumbent_limits(tree)
            return true
        end
    end

    for branch_node in step_obj.branch
        push_to_branch_list!(tree, branch_node)
    end
    return false
end
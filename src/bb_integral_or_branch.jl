"""
    push_integral_or_branch!(m,step_obj,leaf,temp)

Add integral or branch node to step_obj
"""
function push_integral_or_branch!(m,step_obj,leaf,temp)
    # check if all int vars are int
    if are_type_correct(leaf.solution,m.var_type)
        leaf.state = :Integral
        if !temp
            push!(step_obj.integral, leaf)
        end
    else
        leaf.state = :Branch
        if !temp
            push!(step_obj.branch, leaf)
        end
    end
end

"""
    new_integral!(tree,node)

Update the incumbent and add obj constr if in options
Add to solutions if list_of_solutions in options
"""
function new_integral!(tree,node)
    # node.best_bound is the objective for integral values
    tree.nsolutions += 1
    if tree.options.list_of_solutions
        push!(tree.m.solutions, MINLPBnB.SolutionObj(node.solution,node.best_bound))
    end
    if update_incumbent!(tree,node) # returns if new 
        if tree.options.incumbent_constr
            if tree.options.processors > 1
                for p=2:tree.options.processors
                    sendto(p, is_newincumbent=true)
                end
            end
            add_incumbent_constr(tree.m, tree.incumbent)
            tree.m.ncuts += 1
        end
    end
end

"""
    push_to_branch_list!(tree,leaf)

Push a node to the list of branch nodes if better than incumbent
or if all solutions are wanted
"""
function push_to_branch_list!(tree,leaf)
    if tree.options.all_solutions || tree.incumbent == nothing || tree.obj_fac*leaf.best_bound >= tree.obj_fac*tree.incumbent.objval
        push!(tree.branch_nodes, leaf)
    end
end

"""
    upd_integral_branch!(tree,step_obj)

Update the list of integral and branch nodes using the new step_obj
Return true if break
"""
function upd_integral_branch!(tree,step_obj)
    for integral_node in step_obj.integral
        new_integral!(tree,integral_node)
        if isbreak_new_incumbent_limits(tree)
            return true
        end
    end

    for branch_node in step_obj.branch
        push_to_branch_list!(tree,branch_node)
    end
    return false
end
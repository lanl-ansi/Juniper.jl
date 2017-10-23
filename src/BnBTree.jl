module BnBTree
import MINLPBnB
using JuMP



rtol = 1e-6
atol = 1e-6

type BnBNode
    level       :: Int64
    m           :: MINLPBnB.MINLPBnBModel
    children    :: Union{Void,Vector{BnBNode}}
    state       :: Symbol
end

function init(m)
    return BnBNode(1,m,nothing,:Branch)
end

function get_int_variable_idx(num_var,var_type,x)
    println("var_type: ")
    println(var_type)
    println("")
    println("x: ")
    println(x)
    # get variable which is the most unintegral (0.5 isn't integral at all 0.9 is quite a bit :D)
    idx = 1
    max_diff = 0
    for i=1:num_var
        if var_type != :Cont
            diff = abs(x[i]-round(x[i]))
            if diff > max_diff
                idx = i
                max_diff = diff
            end
        end
    end
    println("max_diff: ", max_diff)
    println("idx: ", idx)
    println("x_i: ", x[idx])

    return idx
end

function are_type_correct(sol,types)
    for i=1:length(sol)
        if types[i] != :Cont
            if !isapprox(abs(round(sol[i])-sol[i]),0, atol=atol, rtol=rtol)
                return false
            end
        end
    end
    return true
end

function solve_leaf(leaf)
    status = JuMP.solve(leaf.m.model)
    leaf.m.objval   = getobjectivevalue(leaf.m.model)
    leaf.m.solution = getvalue(leaf.m.x)
    leaf.m.status = status
    if status == :Optimal
        # check if all int vars are int
        if BnBTree.are_type_correct(leaf.m.solution,leaf.m.var_type)
            leaf.state = :Integral
        else
            leaf.state = :Branch
        end
    else
        leaf.state = :Infeasible
    end
    return leaf.state
end

function divide(tree::BnBNode,idx)
    l_m = Base.deepcopy(tree.m)
    r_m = Base.deepcopy(tree.m)

    l_x = l_m.x
    l_cx = l_m.solution[idx]
    r_x = r_m.x
    r_cx = r_m.solution[idx]
    println("Divide")
    @constraint(l_m.model, l_x[idx] <= floor(l_cx))
    @constraint(r_m.model, r_x[idx] >= ceil(r_cx))

    l_tr = BnBNode(tree.level+1,l_m,nothing,:Solve)
    r_tr = BnBNode(tree.level+1,r_m,nothing,:Solve)

    tree.children = [l_tr, r_tr]
    tree.state = :Done

    l_state = solve_leaf(l_tr)
    r_state = solve_leaf(r_tr)
    println("State of left leaf: ", l_state)
    println("State of right leaf: ", r_state)
    println("l sol: ", l_tr.m.solution)
    println("r sol: ", r_tr.m.solution)

    if l_state == :Integral || r_state == :Integral
        # both integral => get better
        if l_state == :Integral && r_state == :Integral
            if l_tr.m.objval > r_tr.m.objval
                return l_tr
            else
                return r_tr
            end
        end
        if l_state == :Integral
            return l_tr
        end
        return r_tr
    end
    # no integral solution ...
    return false    
end

function solve(tree::BnBNode)
    println("Solve Tree")
    # get variable where to split
    m = tree.m
    idx = BnBTree.get_int_variable_idx(m.num_var,m.var_type,m.solution)
    best_leaf = BnBTree.divide(tree,idx)
    return best_leaf.m
end

end
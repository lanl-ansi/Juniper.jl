module BnBTree
import MINLPBnB
using JuMP

rtol = 1e-6
atol = 1e-6

type BnBNode
    parent      :: Union{Void,BnBNode}
    idx         :: Int64
    level       :: Int64
    m           :: MINLPBnB.MINLPBnBModel
    left        :: Union{Void,BnBNode}
    right       :: Union{Void,BnBNode}
    state       :: Symbol
    hasbranchild :: Bool # has child where to branch or is :Branch
end

type BnBTreeObj
    root        :: BnBNode
    incumbent   :: Union{Void,MINLPBnB.MINLPBnBModel}
end

function init(m)
    node = BnBNode(nothing,1,1,m,nothing,nothing,:Branch,true)
    return BnBTreeObj(node,nothing)
end

function get_int_variable_idx(num_var,var_type,x)
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
            leaf.hasbranchild = false
        else
            leaf.state = :Branch
        end
    else
        leaf.state = :Infeasible
        leaf.hasbranchild = false
    end
    return leaf.state
end

function divide!(node::BnBNode,idx)
    l_m = Base.deepcopy(node.m)
    r_m = Base.deepcopy(node.m)

    l_x = l_m.x
    l_cx = l_m.solution[idx]
    r_x = r_m.x
    r_cx = r_m.solution[idx]
    println("Divide")
    @constraint(l_m.model, l_x[idx] <= floor(l_cx))
    @constraint(r_m.model, r_x[idx] >= ceil(r_cx))

    l_tr = BnBNode(node,node.idx*2,node.level+1,l_m,nothing,nothing,:Solve,true)
    r_tr = BnBNode(node,node.idx*2+1,node.level+1,r_m,nothing,nothing,:Solve,true)

    node.left = l_tr
    node.right = r_tr
    node.state = :Done

    l_state = solve_leaf(l_tr)
    r_state = solve_leaf(r_tr)
    println("State of left leaf: ", l_state)
    println("State of right leaf: ", r_state)
    println("l sol: ", l_tr.m.solution)
    println("r sol: ", r_tr.m.solution)

end

function update_incumbent!(tree::BnBTreeObj,node::BnBNode)
    l_nd = node.left
    r_nd = node.right
    l_state, r_state = l_nd.state, r_nd.state
    if l_state == :Integral || r_state == :Integral
        # both integral => get better
        if l_state == :Integral && r_state == :Integral
            if l_nd.m.objval > r_nd.m.objval
                tree.incumbent = l_nd.m
            else
                tree.incumbent = r_nd.m
            end
        elseif l_state == :Integral
            tree.incumbent = l_nd.m
        else
            tree.incumbent = r_nd.m
        end
    end
end

function update_branch!(tree::BnBTreeObj,node::BnBNode)
    l_nd = node.left
    r_nd = node.right
    l_state, r_state = l_nd.state, r_nd.state
    if l_state != :Branch && r_right != :Branch
        node.hasbranchild = false
   
        # both children aren't branch nodes
        # bubble up to check where to set node.hasbranchild = false
        while node.parent != nothing
            node = node.parent
            if node.left.hasbranchild || node.right.hasbranchild
                break
            else
                node.hasbranchild = false
            end
        end
    end
end

function get_best_branch_node(tree::BnBTreeObj)
    node = tree.root
    obj_sense = tree.root.m.obj_sense
    factor = 1
    if obj_sense == :Min
        factor = -1
    end

    if node.state == :Branch
        return node
    end

    while true
        l_nd = node.left
        r_nd = node.right
        if node.hasbranchild == true
            if l_nd.state == :Branch && r_nd.state == :Branch 
                # use node with better obj
                if factor*l_nd.m.objval > factor*r_nd.m.objval
                    return l_nd
                else
                    return r_nd
                end
            elseif l_nd.state == :Branch
                return l_nd
            elseif r_nd.state == :Branch
                return r_nd
            else
                # get into better branch
                if l_nd.hasbranchild && r_nd.hasbranchild
                    exit("get into better branch")
                elseif !l_nd.hasbranchild && !r_nd.hasbranchild
                    exit("Infeasible")
                elseif l_nd.hasbranchild
                    node = l_nd
                else
                    node = r_nd
                end
            end
        end
    end
end

function print(node::BnBNode)
    indent = (node.level-1)*2
    indent_str = ""
    for i=1:indent
        indent_str *= " "
    end
    println(indent_str*"idx"*": "*string(node.idx))
    println(indent_str*"state"*": "*string(node.state))
    println(indent_str*"hasbranchild"*": "*string(node.hasbranchild))
end

function print_rec(node::BnBNode)
    print(node)
    if node.left != nothing
        print_rec(node.left)
    end
    if node.right != nothing
        print_rec(node.right)
    end
end

function print(tree::BnBTreeObj)
    node = tree.root
    print_rec(node)
end

function solve(tree::BnBTreeObj)
    println("Solve Tree")
    # get variable where to split
    node = tree.root

    while true
        println("=================================")
        println("NodeIdx: ", node.idx)
        println("=================================")
        m = node.m
        v_idx = BnBTree.get_int_variable_idx(m.num_var,m.var_type,m.solution)
        println("v_idx: ", v_idx)

        BnBTree.divide!(node,v_idx)
        print(tree)

        # update incumbent
        BnBTree.update_incumbent!(tree,node)
        if tree.incumbent != nothing
            return tree.incumbent
        end

        BnBTree.update_branch!(tree,node)
        # get best branch node
        node = BnBTree.get_best_branch_node(tree)
        println("next nodeidx: ", node.idx)
        
        if node.level == 4
            break
        end
    end
end

end
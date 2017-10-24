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
    best_bound :: Union{Void,Float64}
end

type BnBTreeObj
    root        :: BnBNode
    incumbent   :: Union{Void,MINLPBnB.MINLPBnBModel}
    print_syms  :: Vector{Symbol}
end

function init(m)
    node = BnBNode(nothing,1,1,m,nothing,nothing,:Branch,true,m.objval)
    return BnBTreeObj(node,nothing,m.print_syms)
end

function check_print(vec::Vector{Symbol}, ps::Vector{Symbol})
    for v in vec
        if v in ps
            return true
        end
    end
    return false
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
            leaf.best_bound = leaf.m.objval
        else
            leaf.state = :Branch
            leaf.best_bound = leaf.m.objval
        end
    else
        leaf.state = :Infeasible
        leaf.hasbranchild = false
    end
    return leaf.state
end

function branch!(node::BnBNode,idx,ps)
    l_m = Base.deepcopy(node.m)
    r_m = Base.deepcopy(node.m)

    l_x = l_m.x
    l_cx = l_m.solution[idx]
    r_x = r_m.x
    r_cx = r_m.solution[idx]
    BnBTree.check_print(ps,[:All,:FuncCall]) && println("branch")
    @constraint(l_m.model, l_x[idx] <= floor(l_cx))
    @constraint(r_m.model, r_x[idx] >= ceil(r_cx))

    l_tr = BnBNode(node,node.idx*2,node.level+1,l_m,nothing,nothing,:Solve,true,nothing)
    r_tr = BnBNode(node,node.idx*2+1,node.level+1,r_m,nothing,nothing,:Solve,true,nothing)

    node.left = l_tr
    node.right = r_tr
    node.state = :Done

    l_state = solve_leaf(l_tr)
    r_state = solve_leaf(r_tr)

    if BnBTree.check_print(ps,[:All])
        println("State of left leaf: ", l_state)
        println("State of right leaf: ", r_state)
        println("l sol: ", l_tr.m.solution)
        println("r sol: ", r_tr.m.solution)
    end

end

function update_incumbent!(tree::BnBTreeObj,node::BnBNode)
    ps = tree.print_syms
    BnBTree.check_print(ps,[:All,:FuncCall]) && println("update_incumbent")

    l_nd = node.left
    r_nd = node.right
    l_state, r_state = l_nd.state, r_nd.state
    factor = 1
    if tree.root.m.obj_sense == :Min
        factor = -1
    end

    if l_state == :Integral || r_state == :Integral
        # both integral => get better
        if l_state == :Integral && r_state == :Integral
            if factor*l_nd.m.objval > factor*r_nd.m.objval
                possible_incumbent = l_nd.m
            else
                possible_incumbent = r_nd.m
            end
        elseif l_state == :Integral
            possible_incumbent = l_nd.m
        else
            possible_incumbent = r_nd.m
        end
        if tree.incumbent == nothing || factor*possible_incumbent.objval > factor*tree.incumbent.objval
            tree.incumbent = possible_incumbent
            return true
        end
    end
 
    return false
end

function update_branch!(tree::BnBTreeObj,node::BnBNode)
    ps = tree.print_syms
    l_nd = node.left
    r_nd = node.right
    l_state, r_state = l_nd.state, r_nd.state
    BnBTree.check_print(ps,[:All,:FuncCall]) && println("update branch")
    BnBTree.check_print(ps,[:All]) && println(l_state, " ", r_state)
    if l_state != :Branch && r_state != :Branch
        node.hasbranchild = false
   
        # both children aren't branch nodes
        # bubble up to check where to set node.hasbranchild = false
        
        while node.parent != nothing
            node = node.parent
            BnBTree.check_print(ps,[:All]) && println("node.level: ", node.level)
            if node.left.hasbranchild || node.right.hasbranchild
                BnBTree.check_print(ps,[:All]) && println("break")
                break
            else
                node.hasbranchild = false
            end
        end
    else
        factor = 1
        if tree.root.m.obj_sense == :Min
            factor = -1
        end

        while node != nothing
            BnBTree.check_print(ps,[:All]) && println("Node idx: ", node.idx)
            l_nd = node.left
            r_nd = node.right
            
            if l_nd.best_bound == nothing
                node.best_bound = r_nd.best_bound
            elseif r_nd.best_bound == nothing
                node.best_bound = l_nd.best_bound
            elseif factor*l_nd.best_bound > factor*r_nd.best_bound
                node.best_bound = l_nd.best_bound
            else
                node.best_bound = r_nd.best_bound
            end
            node = node.parent
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
                    if factor*l_nd.best_bound > factor*r_nd.best_bound
                        node = l_nd
                    else
                        node = r_nd
                    end
                elseif !l_nd.hasbranchild && !r_nd.hasbranchild
                    println("node idx: ", node.idx)
                    print(tree;remove=:hasnobranchild)
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

function prune!(node::BnBNode, value)
    obj_sense = node.m.obj_sense
    factor = 1
    if obj_sense == :Min
        factor = -1
    end
    if node.hasbranchild && factor*value > factor*node.best_bound
        node.hasbranchild = false 
    else
        if node.left != nothing
            prune!(node.left, value)
        end
        if node.right != nothing
            prune!(node.left, value)
        end
    end
end

function prune!(tree::BnBTreeObj)
    incumbent_val = tree.incumbent.objval
    ps = tree.print_syms
    BnBTree.check_print(ps,[:All,:Incumbent]) && println("incumbent_val: ", incumbent_val)

    prune!(tree.root, incumbent_val)
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
    println(indent_str*"best_bound"*": "*string(node.best_bound))
end

function print_rec(node::BnBNode;remove=false)
    if remove != :hasnobranchild || node.hasbranchild
        print(node)
        if node.left != nothing
            print_rec(node.left;remove=remove)
        end
        if node.right != nothing
            print_rec(node.right;remove=remove)
        end
    end
end

function print(tree::BnBTreeObj;remove=false)
    node = tree.root
    print_rec(node;remove=remove)
end

function solve(tree::BnBTreeObj)
    ps = tree.print_syms
    BnBTree.check_print(ps,[:All,:FuncCall]) && println("Solve Tree")
    # get variable where to split
    node = tree.root
    

    while true
        m = node.m
        v_idx = BnBTree.get_int_variable_idx(m.num_var,m.var_type,m.solution)
        BnBTree.check_print(ps,[:All]) && println("v_idx: ", v_idx)

        BnBTree.branch!(node,v_idx,tree.print_syms)
    
        BnBTree.update_branch!(tree,node)

        # update incumbent
        if BnBTree.update_incumbent!(tree,node)
            BnBTree.check_print(ps,[:All]) && println("Prune")
            BnBTree.prune!(tree)
            BnBTree.check_print(ps,[:All]) && println("pruned")            
            BnBTree.check_print(ps,[:All]) && print(tree)
            BnBTree.check_print(ps,[:All,:NewIncumbent]) && println("New incumbent: ",tree.incumbent.objval)
        end
        # check if best
        if !tree.root.hasbranchild
            return tree.incumbent
        end
        # println("Best bound: ", tree.root.best_bound)
        # println("Node level: ", node.level)

        # print(tree)
        # get best branch node
        node = BnBTree.get_best_branch_node(tree)
        
    end
end

end
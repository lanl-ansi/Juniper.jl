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
    var_idx     :: Int64
    left        :: Union{Void,BnBNode}
    right       :: Union{Void,BnBNode}
    state       :: Symbol
    hasbranchild :: Bool # has child where to branch or is :Branch
    best_bound  :: Union{Void,Float64}
end

type BnBTreeObj
    root        :: BnBNode
    incumbent   :: Union{Void,MINLPBnB.MINLPBnBModel}
    print_syms  :: Vector{Symbol}
    obj_gain    :: Vector{Float64} # gain of objective per variable
    obj_gain_c  :: Vector{Float64} # obj_gain / obj_gain_c => average gain
    int2var_idx :: Vector{Int64}
    var2int_idx :: Vector{Int64}
end

function init(m)
    node = BnBNode(nothing,1,1,m,0,nothing,nothing,:Branch,true,m.objval)
    obj_gain = zeros(m.num_int_bin_var)
    obj_gain_c = ones(m.num_int_bin_var)
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
    return BnBTreeObj(node,nothing,m.print_syms,obj_gain,obj_gain_c,int2var_idx,var2int_idx)
end

function new_default_node(parent,idx,level,m;
                            var_idx=0,left=nothing,right=nothing,
                            state=:Solve,hasbranchild=true,best_bound=nothing)

    return BnBNode(parent,idx,level,m,var_idx,left,right,state,hasbranchild,best_bound)     
end

function check_print(vec::Vector{Symbol}, ps::Vector{Symbol})
    for v in vec
        if v in ps
            return true
        end
    end
    return false
end

function branch_mostinfeasible(tree,node,num_var,var_type,x)
    # get variable which is the most unintegral (0.5 isn't integral at all 0.9 is quite a bit :D)
    idx = 0
    max_diff = 0
    for i=1:num_var
        if var_type[i] != :Cont
            diff = abs(x[i]-round(x[i]))
            if diff > max_diff
                idx = i
                max_diff = diff
            end
        end
    end
    return idx
end

"""
    get_int_variable_idx(tree,node,branch_strat,num_var,var_type,x;counter=1)

Get the index of a variable to branch on.
branch_strat might be :MostInfeasible or :PseudoCost
"""
function get_int_variable_idx(tree,node,branch_strat,num_var,var_type,x;counter=1)
    # Most Infeasible Branching
    idx = 0
    if branch_strat == :MostInfeasible
        return BnBTree.branch_mostinfeasible(tree,node,num_var,var_type,x)
    elseif branch_strat == :PseudoCost
        if counter == 1
            idx = BnBTree.branch_mostinfeasible(tree,node,num_var,var_type,x)
        else
            # use the one with highest obj_gain which is currently continous
            obj_gain_average = tree.obj_gain./tree.obj_gain_c
            sort_idx = tree.int2var_idx[sortperm(obj_gain_average, rev=true)]
            for l_idx in sort_idx
                if !is_type_correct(x[l_idx],var_type[l_idx])
                    sol = node.m.solution[l_idx]
                    u_b = node.m.u_var[l_idx]
                    l_b = node.m.l_var[l_idx]
                    if isapprox(u_b,floor(sol),atol=atol) || isapprox(l_b, ceil(sol),atol=atol)
                        continue
                    end
                    return l_idx
                end
            end
        end
    end
    @assert idx != 0
    return idx
end

"""
    is_type_correct(x,var_type)

Check whether a variable x has the correct type
"""
function is_type_correct(x,var_type)
    if var_type != :Cont
        if !isapprox(abs(round(x)-x),0, atol=atol, rtol=rtol)
           return false
        end
    end
    return true
end

"""
    are_type_correct(sol,types)

Check whether all variables have the correct type
"""
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

"""
    solve_leaf(leaf)

Solve a leaf by relaxation leaf is just a node.
Set the state,hasbranchild and best_bound property
Return state
"""
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

"""
    branch!(node::BnBNode,idx,ps)

Branch a node by using x[idx] <= floor(x[idx]) and x[idx] >= ceil(x[idx])
Solve both nodes and set current node state to done.
"""
function branch!(node::BnBNode,idx,ps)
    l_m = Base.deepcopy(node.m)
    r_m = Base.deepcopy(node.m)

    # save that this node branches on this particular variable
    node.var_idx = idx

    l_x = l_m.x
    l_cx = l_m.solution[idx]
    r_x = r_m.x
    r_cx = r_m.solution[idx]
    BnBTree.check_print(ps,[:All,:FuncCall]) && println("branch")
    
    if isapprox(l_m.u_var[idx],floor(l_cx),atol=atol) || isapprox(r_m.l_var[idx], ceil(r_cx),atol=atol)
        error("Shouldn't solve again")
    end
    @constraint(l_m.model, l_x[idx] <= floor(l_cx))    
    l_m.u_var[idx] = floor(l_cx)
    @constraint(r_m.model, r_x[idx] >= ceil(r_cx))
    r_m.l_var[idx] = ceil(r_cx)

    l_tr = BnBTree.new_default_node(node,node.idx*2,node.level+1,l_m)
    r_tr = BnBTree.new_default_node(node,node.idx*2+1,node.level+1,r_m)

    node.left = l_tr
    node.right = r_tr
    node.state = :Done

    leaf_start = time()
    l_state = solve_leaf(l_tr)
    r_state = solve_leaf(r_tr)
    leaf_time = time()-leaf_start

    if BnBTree.check_print(ps,[:All])
        println("State of left leaf: ", l_state)
        println("State of right leaf: ", r_state)
        println("l sol: ", l_tr.m.solution)
        println("r sol: ", r_tr.m.solution)
    end
    return leaf_time
end

"""
    update_gains(tree::BnBTreeObj,node::BnBNode;counter=1)

Update the objective gains for the branch variable used for node
"""
function update_gains(tree::BnBTreeObj,node::BnBNode;counter=1)
    gain = 0
    frac_val = node.m.solution[node.var_idx]
    if node.left.state == :Branch
        int_val = floor(frac_val)
        gain += abs(node.best_bound-node.left.best_bound)/abs(frac_val-int_val)
    end
    if node.right.state == :Branch
        int_val = ceil(frac_val)
        gain += abs(node.best_bound-node.right.best_bound)/abs(frac_val-int_val)
    end

    # update all (just average of the one branch we have)
    if counter == 1
        tree.obj_gain += gain
    else
        idx = tree.var2int_idx[node.var_idx]
        tree.obj_gain[idx] += gain
        tree.obj_gain_c[idx] += 1
    end
end

"""
    update_incumbent!(tree::BnBTreeObj,node::BnBNode)

Update the incumbent if there is a new Integral solution which is better.
Return true if updated false otherwise
"""
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

"""
    update_branch!(tree::BnBTreeObj,node::BnBNode)

Update the branch tree. If on both children can't be branched on
=> set hasbranchild = false and check the parents as well (bubble up)
If one of both children can be branched on => bubble up the best_bound
"""
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
        # Bubble up the best bound of the children
        # => The root has always the best bound of all of it's children
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

"""
    get_best_branch_node(tree::BnBTreeObj)

Get the index of the breach node which should be used for the next branch.
Currently get's the branch with the worst best bound
"""
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
                # use node with worse obj
                if factor*l_nd.m.objval < factor*r_nd.m.objval
                    return l_nd
                else
                    return r_nd
                end
            elseif l_nd.state == :Branch
                return l_nd
            elseif r_nd.state == :Branch
                return r_nd
            else
                # get into worse branch
                if l_nd.hasbranchild && r_nd.hasbranchild
                    if factor*l_nd.best_bound < factor*r_nd.best_bound
                        node = l_nd
                    else
                        node = r_nd
                    end
                elseif !l_nd.hasbranchild && !r_nd.hasbranchild
                    println("node idx: ", node.idx)
                    print(tree)
                    error("Infeasible")
                elseif l_nd.hasbranchild
                    node = l_nd
                else
                    node = r_nd
                end
            end
        end
    end
end

"""
    prune!(node::BnBNode, value)

Get rid of nodes which have a worse best bound then specified by value. 
Is recursive
"""
function prune!(node::BnBNode, value)
    obj_sense = node.m.obj_sense
    factor = 1
    if obj_sense == :Min
        factor = -1
    end
    if node.hasbranchild && factor*value > factor*node.best_bound
        node.hasbranchild = false 
        node.left = nothing
        node.right = nothing
    else
        if node.left != nothing
            prune!(node.left, value)
        end
        if node.right != nothing
            prune!(node.left, value)
        end
    end
end

"""
    prune!(tree::BnBTreeObj)

Call prune! for the root node using the incumbent value
"""
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
    println(indent_str*"var_idx"*": "*string(node.var_idx))
    println(indent_str*"state"*": "*string(node.state))
    # println(indent_str*"hasbranchild"*": "*string(node.hasbranchild))
    # println(indent_str*"best_bound"*": "*string(node.best_bound))
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

"""
    solve(tree::BnBTreeObj)

Solve the MIP part of a problem given by BnBTreeObj using branch and bound.
 - Identify the node to branch on
 - Get variable to branch on
 - Solve subproblems
"""
function solve(tree::BnBTreeObj)
    ps = tree.print_syms
    BnBTree.check_print(ps,[:All,:FuncCall]) && println("Solve Tree")
    # get variable where to split
    node = tree.root
    counter = 1    

    branch_strat = :PseudoCost
    time_upd_gains = 0
    time_get_idx = 0
    time_branch = 0
    time_solve_leafs = 0
    
    time_bnb_solve_start = time()
    while true
        println("Best bound: ", tree.root.best_bound)
        m = node.m
        println("Node level: ", node.level)
        get_idx_start = time()
        v_idx = BnBTree.get_int_variable_idx(tree,node,branch_strat, m.num_var,m.var_type,m.solution;counter=counter)
        time_get_idx += time()-get_idx_start
    
        BnBTree.check_print(ps,[:All]) && println("v_idx: ", v_idx)

        branch_start = time()
        time_solve_leafs += BnBTree.branch!(node,v_idx,tree.print_syms)
        time_branch += time()-branch_start

        if branch_strat == :PseudoCost
            upd_start = time()
            BnBTree.update_gains(tree,node;counter=counter)    
            time_upd_gains += time()-upd_start
        end

        BnBTree.update_branch!(tree,node)

        # update incumbent
        if BnBTree.update_incumbent!(tree,node)
            BnBTree.check_print(ps,[:All]) && println("Prune")
            BnBTree.prune!(tree)
            BnBTree.check_print(ps,[:All]) && println("pruned")            
            BnBTree.check_print(ps,[:All]) && print(tree)
            if BnBTree.check_print(ps,[:All,:NewIncumbent]) 
                println("New incumbent: ",tree.incumbent.objval)
                println("Best bound: ",tree.root.best_bound)
            end
        end
        # check if best
        if !tree.root.hasbranchild
            break
        end
        # println("Best bound: ", tree.root.best_bound)
        # println("Node level: ", node.level)

        # print(tree)
        # get best branch node
        node = BnBTree.get_best_branch_node(tree)

        counter += 1
    end

    # print(tree)
   
    time_bnb_solve = time()-time_bnb_solve_start
    println("#branches: ", counter)
    println("BnB time: ", round(time_bnb_solve,2))
    println("Solve leaf time: ", round(time_solve_leafs,2))
    println("Branch time: ", round(time_branch,2))
    println("Get idx time: ", round(time_get_idx,2))
    println("Upd gains time: ", round(time_upd_gains,2))
    return tree.incumbent
end

end
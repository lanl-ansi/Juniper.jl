module BnBTree
import MINLPBnB
using JuMP
using Ipopt
using MathProgBase

rtol = 1e-6
atol = 1e-6
time_solve_leaf!s_get_idx = 0.0
time_solve_leaf!s_branch = 0.0

type BnBNode
    idx         :: Int64
    level       :: Int64
    l_var       :: Vector{Float64}
    u_var       :: Vector{Float64}
    solution    :: Vector{Float64}
    var_idx     :: Int64
    state       :: Symbol
    best_bound  :: Union{Void,Float64}
end

type IncumbentSolution
    objval      :: Float64
    solution    :: Vector{Float64}
    status      :: Symbol
    best_bound  :: Float64
end

type BnBTreeObj
    m           :: MINLPBnB.MINLPBnBModel
    incumbent   :: Union{Void,IncumbentSolution}
    obj_gain    :: Vector{Float64} # gain of objective per variable
    obj_gain_c  :: Vector{Float64} # obj_gain / obj_gain_c => average gain
    int2var_idx :: Vector{Int64}
    var2int_idx :: Vector{Int64}
    options     :: MINLPBnB.SolverOptions
    obj_fac     :: Int64 # factor for objective 1 if max -1 if min
    start_time  :: Float64 
    nsolutions  :: Int64
    branch_nodes:: Vector{Union{Void,BnBNode}}
    best_bound  :: Float64
end

# the object holds information for the current step
type StepObj
    node        :: BnBNode # current branch node
    var_idx     :: Int64   # variable to branch on
    state       :: Symbol  # if infeasible => break (might be set by strong branching)
    nrestarts   :: Int64 
    idx_time    :: Float64
    leaf_idx_time :: Float64
    upd_gains_time :: Float64
    leaf_branch_time :: Float64
    branch_time     :: Float64
end

function init(start_time, m)
    srand(1)
    node = BnBNode(1,1,m.l_var,m.u_var,m.solution,0,:Branch,m.objval)
    obj_gain = zeros(m.num_int_bin_var)
    obj_gain_c = zeros(m.num_int_bin_var)
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
    factor = 1
    if m.obj_sense == :Min
        factor = -1
    end
    return BnBTreeObj(m,nothing,obj_gain,obj_gain_c,int2var_idx,var2int_idx,m.options,factor,start_time,0,[node],NaN)
end

function new_default_node(idx,level,l_var,u_var,solution;
                            var_idx=0,
                            state=:Solve,best_bound=nothing)

    l_var = copy(l_var)
    u_var = copy(u_var)
    solution = copy(solution)
    return BnBNode(idx,level,l_var,u_var,solution,var_idx,state,best_bound)     
end

function new_default_step_obj(node)
    idx_time = 0.0
    leaf_idx_time = 0.0
    upd_gains_time = 0.0
    leaf_branch_time = 0.0
    branch_time = 0.0
    return StepObj(node,0,:None,0,idx_time,leaf_idx_time,upd_gains_time,leaf_branch_time,branch_time)
end

function check_print(vec::Vector{Symbol}, ps::Vector{Symbol})
    for v in vec
        if v in ps
            return true
        end
    end
    return false
end

"""
    branch_mostinfeasible(tree,node)

Get the index of an integer variable which is currently continuous which is most unintegral.
(nearest to *.5)
"""
function branch_mostinfeasible(tree,node)
    x = node.solution
    idx = 0
    max_diff = 0
    for i=1:tree.m.num_var
        if tree.m.var_type[i] != :Cont
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
    init_strong_restart!(node, var_idx, int_var_idx, l_nd, r_nd, reasonable_int_vars, infeasible_int_vars)

Tighten the bounds for the node and check if there are variables that need to be checked for a restart.
"""
function init_strong_restart!(node, var_idx, int_var_idx, l_nd, r_nd, reasonable_int_vars, infeasible_int_vars)
    restart = false

    # set the bounds directly for the node
    # also update the best bound and the solution
    if l_nd.state == :Infeasible
        node.l_var[var_idx] = ceil(node.solution[var_idx])
        node.best_bound = r_nd.best_bound
        node.solution = r_nd.solution
    else
        node.u_var[var_idx] = floor(node.solution[var_idx])
        node.best_bound = l_nd.best_bound
        node.solution = l_nd.solution
    end
    
    push!(infeasible_int_vars,int_var_idx)

    if length(reasonable_int_vars) == length(infeasible_int_vars)
        # basically branching on the last infeasible variable 
        max_gain_var = infeasible_int_vars[end]
        strong_int_vars = [infeasible_int_vars[end]] # don't divide by 0 later
        node.left = l_nd
        node.right = r_nd
    else
        max_gain_var = 0
        strong_int_vars = zeros(Int64,0)
        restart = true
    end
    return restart, infeasible_int_vars, max_gain_var, strong_int_vars
end

"""
    branch_strong(tree,node,counter)

Try to branch on a few different variables and choose the one with highest obj_gain.
Update obj_gain for the variables tried and average the other ones.
"""
function branch_strong(tree,step_obj,counter)
    function init_variables()
        max_gain = 0.0
        max_gain_var = 0
        strong_int_vars = zeros(Int64,0)
        return max_gain, max_gain_var, strong_int_vars
    end
    node = step_obj.node

    strong_restarts = -1 

    # generate an of variables to branch on
    num_strong_var = tree.options.strong_branching_nvars

    # get reasonable candidates (not type correct and not already perfectly bounded)
    int_vars = tree.m.num_int_bin_var
    reasonable_int_vars = zeros(Int64,0)
    for i=1:int_vars
        idx = tree.int2var_idx[i]
        u_b = node.u_var[idx]
        l_b = node.l_var[idx]
        if isapprox(u_b,l_b,atol=atol) || BnBTree.is_type_correct(node.solution[idx],tree.m.var_type[idx])
            continue
        end
        push!(reasonable_int_vars,i)
    end
    shuffle!(reasonable_int_vars)
    reasonable_int_vars = reasonable_int_vars[1:minimum([num_strong_var,length(reasonable_int_vars)])]

    # compute the gain for each reasonable candidate and choose the highest
    max_gain, max_gain_var, strong_int_vars = init_variables()
    left_node = nothing
    right_node = nothing

    restart = true
    infeasible_int_vars = zeros(Int64,0)
    status = :Normal
    while restart 
        strong_restarts += 1
        restart = false
        for int_var_idx in reasonable_int_vars
            # don't rerun if the variable has already one infeasible node
            if int_var_idx in infeasible_int_vars
                continue
            end
            push!(strong_int_vars, int_var_idx)
            var_idx = tree.int2var_idx[int_var_idx]
            step_obj.var_idx = var_idx
            u_b, l_b = node.u_var[var_idx], node.l_var[var_idx]
            # don't rerun if bounds are exact or is type correct
            if isapprox(u_b,l_b,atol=atol) || BnBTree.is_type_correct(node.solution[var_idx],tree.m.var_type[var_idx])
                continue
            end
            # branch on the current variable and get the corresponding children
            l_nd,r_nd = BnBTree.branch!(tree,step_obj,counter;temp=true)

            if l_nd.state == :Infeasible && r_nd.state == :Infeasible && counter == 1
                status = :Infeasible
                break
            end

            # if restart is true => check if one part is infeasible => update bounds & restart
            if tree.options.strong_restart == true
                if l_nd.state == :Infeasible || r_nd.state == :Infeasible
                    max_gain = 0.0
                    if l_nd.state == :Infeasible && r_nd.state == :Infeasible
                        node.state = :Infeasible
                        break
                    end
                    restart,infeasible_int_vars,max_gain_var,strong_int_vars = BnBTree.init_strong_restart!(node, var_idx, int_var_idx, l_nd, r_nd, reasonable_int_vars, infeasible_int_vars)
                end
            end
            gain = BnBTree.compute_gain(node;l_nd=l_nd,r_nd=r_nd,inf=true)
            if gain > max_gain
                max_gain = gain
                max_gain_var = var_idx
                left_node = l_nd
                right_node = r_nd
                # gain is set to inf if Integral or Infeasible
                # TODO: Might be reasonable to use something different
                if gain == Inf
                    break
                end
            end
            tree.obj_gain[int_var_idx] = gain
        end
    end
    
    if status != :Infeasible
            for leaf in [left_node,right_node]
                if leaf.state == :Branch
                    BnBTree.push_to_branch_list!(tree,leaf)
                elseif leaf.state == :Integral
                    BnBTree.new_integral!(tree,leaf)
                end
            end   

        # set the variable to branch (best gain)
        node.state = :Done
        node.var_idx = max_gain_var

        # assign values to untested variables only for the first strong branch
        if counter == 1
            # all other variables that haven't been checked get the median value of the others
            med_gain = median(tree.obj_gain[strong_int_vars])
            rest = filter(i->!(i in strong_int_vars),1:int_vars)
            tree.obj_gain[rest] += med_gain
            tree.obj_gain_c += 1
        else
            tree.obj_gain_c[strong_int_vars] += 1
        end
    end

    @assert max_gain_var != 0 || status == :Infeasible || node.state == :Infeasible
    return status, max_gain_var, strong_restarts
end

"""
    upd_int_variable_idx!(step_obj,tree,node,counter=1)    

Get the index of a variable to branch on.
"""
function upd_int_variable_idx!(step_obj,tree,counter::Int64=1)   
    start = time()
    node = step_obj.node
    idx = 0
    strong_restarts = 0
    branch_strat = tree.options.branch_strategy
    status = :Normal
    if branch_strat == :MostInfeasible
        idx = BnBTree.branch_mostinfeasible(tree,node)
    elseif branch_strat == :PseudoCost || branch_strat == :StrongPseudoCost
        if counter == 1 && branch_strat == :PseudoCost
            idx = BnBTree.branch_mostinfeasible(tree,node)
        elseif counter <= tree.options.strong_branching_nsteps && branch_strat == :StrongPseudoCost
            status, idx, strong_restarts = BnBTree.branch_strong(tree,step_obj,counter)
        else
            # use the one with highest obj_gain which is currently continous
            obj_gain_average = tree.obj_gain./tree.obj_gain_c
            sort_idx = tree.int2var_idx[sortperm(obj_gain_average, rev=true)]
            for l_idx in sort_idx
                if !is_type_correct(node.solution[l_idx],tree.m.var_type[l_idx])
                    u_b = node.u_var[l_idx]
                    l_b = node.l_var[l_idx]
                    # if the upper bound is the lower bound => no reason to branch
                    if isapprox(u_b,l_b,atol=atol)
                        continue
                    end
                    idx = l_idx
                    break
                end
            end
        end
    end
    step_obj.state = status
    step_obj.var_idx = idx
    step_obj.nrestarts = strong_restarts
    step_obj.idx_time = time()-start
    return
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
    solve_leaf!(tree,leaf)

Solve a leaf by relaxation leaf is just a node.
Set the state and best_bound property
Update incumbent if new and add node to branch list if :Branch
Return state
"""
function solve_leaf!(tree,leaf,temp)
     # set bounds
    for i=1:tree.m.num_var
        JuMP.setlowerbound(tree.m.x[i], leaf.l_var[i])    
        JuMP.setupperbound(tree.m.x[i], leaf.u_var[i])
    end

    status = JuMP.solve(tree.m.model)
    objval   = getobjectivevalue(tree.m.model)
    leaf.solution = getvalue(tree.m.x)
    status = status
    if status == :Error
        # println(leaf.m.model)
        println(Ipopt.ApplicationReturnStatus[internalmodel(tree.m.model).inner.status])
        # error("...")
        leaf.state = :Error
    elseif status == :Optimal
        # check if all int vars are int
        if BnBTree.are_type_correct(leaf.solution,tree.m.var_type)
            leaf.state = :Integral
            leaf.best_bound = objval
            if !temp
                if BnBTree.new_integral!(tree,leaf) == :Break
                    leaf.state = :Break
                end
            end
        else
            leaf.state = :Branch
            leaf.best_bound = objval
            if !temp
                push_to_branch_list!(tree,leaf)
            end
        end
    else
        leaf.state = :Infeasible
    end
    return leaf.state
end

function new_integral!(tree,node)
    # node.best_bound is the objective for integral values
    tree.nsolutions += 1
    if tree.options.list_of_solutions
        push!(tree.m.solutions, MINLPBnB.SolutionObj(node.solution,node.best_bound))
    end
    if update_incumbent!(tree,node) # returns if new 
        BnBTree.add_obj_constr(tree)
    end
    if break_new_incumbent_limits(tree)
        return :Break
    end
end

function push_to_branch_list!(tree,leaf)
    if tree.incumbent == nothing || tree.obj_fac*leaf.best_bound >= tree.obj_fac*tree.incumbent.objval
        # tree.incumbent != nothing && println("Inc", tree.incumbent.objval)
        # println("Push to list: ", leaf.best_bound)
        push!(tree.branch_nodes, leaf)
    end
end

"""
    branch!(tree,step_obj,counter)

Branch a node by using x[idx] <= floor(x[idx]) and x[idx] >= ceil(x[idx])
Solve both nodes and set current node state to done.
"""
function branch!(tree,step_obj,counter;temp=false)
    global time_solve_leaf!s_get_idx, time_solve_leaf!s_branch
    ps = tree.options.log_levels
    node = step_obj.node
    vidx = step_obj.var_idx

    start = time()

    # it might be already branched on
    if node.state != :Branch
        return nothing,nothing
    end
    
    l_nd = BnBTree.new_default_node(node.idx*2,node.level+1,node.l_var,node.u_var,node.solution)
    r_nd = BnBTree.new_default_node(node.idx*2+1,node.level+1,node.l_var,node.u_var,node.solution)

    l_nd.u_var[vidx] = floor(node.solution[vidx])
    r_nd.l_var[vidx] = ceil(node.solution[vidx])

    # save that this node branches on this particular variable
    node.var_idx = vidx

    BnBTree.check_print(ps,[:All,:FuncCall]) && println("branch")
    
    if !temp
        node.state = :Done
    end
    
    start_leaf = time()
    l_state = solve_leaf!(tree,l_nd,temp)
    r_state = solve_leaf!(tree,r_nd,temp)
    leaf_time = time() - start_leaf

    if temp
        step_obj.leaf_idx_time += leaf_time
    else
        step_obj.leaf_branch_time += leaf_time
    end

    if l_state == :Break || r_state == :Break
        step_obj.state = :Break
    end

    branch_strat = tree.options.branch_strategy
    if branch_strat == :PseudoCost || (branch_strat == :StrongPseudoCost && counter > tree.options.strong_branching_nsteps)
        upd_start = time()
        gain_gap = BnBTree.update_gains!(tree,node,l_nd,r_nd,counter)    
        step_obj.upd_gains_time = time()-upd_start
    end

    if BnBTree.check_print(ps,[:All])
        println("State of left leaf: ", l_state)
        println("State of right leaf: ", r_state)
        println("l sol: ", l_nd.solution)
        println("r sol: ", r_nd.solution)
    end

    if !temp
        step_obj.branch_time += time()-start
    end

    return l_nd,r_nd
end

"""
    compute_gain(node;l_nd::BnBNode=node.left,r_nd::BnBNode=node.right)

Compute the gain of the children to it's parent.
gain = abs(ob-nb)/abs(frac_val-int_val) where ob = old best bound, nb = new best bound
If the state of one child is Integral or Infeasible
    return Inf
else return the smaller gain of both children
"""
function compute_gain(node;l_nd::BnBNode=node.left,r_nd::BnBNode=node.right,inf=false)
    gc = 0
    gain_l = 0.0
    gain_r = 0.0
    frac_val = node.solution[node.var_idx]
    if inf && (l_nd.state == :Integral || r_nd.state == :Integral || l_nd.state == :Infeasible || r_nd.state == :Infeasible)
        return Inf
    end
    if l_nd.state == :Error && r_nd.state == :Error
        return 0.0
    end
    if l_nd.state == :Branch || l_nd.state == :Integral
        int_val = floor(frac_val)
        gain_l = abs(node.best_bound-l_nd.best_bound)/abs(frac_val-int_val) 
        gc += 1
    end
    if r_nd.state == :Branch || r_nd.state == :Integral
        int_val = ceil(frac_val)
        gain_r = abs(node.best_bound-r_nd.best_bound)/abs(frac_val-int_val)
        gc += 1
    end
    gc == 0 && return 0.0
    # use always minimum of both
    return gain_r > gain_l ? gain_l : gain_r
end

"""
    update_gains!(tree::BnBTreeObj,node::BnBNode,counter)

Update the objective gains for the branch variable used for node
"""
function update_gains!(tree::BnBTreeObj,parent::BnBNode,l_nd,r_nd,counter)
    gain = BnBTree.compute_gain(parent;l_nd=l_nd,r_nd=r_nd)

    idx = tree.var2int_idx[parent.var_idx]
    guess = tree.obj_gain[idx]/tree.obj_gain_c[idx]
    if gain == 0 && guess == 0
        gap = 0.0
    elseif gain == 0
        gap = Inf
    else
        gap = abs(guess-gain)/gain*100    
    end

    # update all (just average of the one branch we have)
    if counter == 1
        tree.obj_gain += gain
    else
        tree.obj_gain[idx] += gain
        tree.obj_gain_c[idx] += 1
    end
    return gap
end

"""
    update_incumbent!(tree::BnBTreeObj,node::BnBNode)

Get's called if new integral solution was found. 
Check whether it's a new incumbent and update if necessary
"""
function update_incumbent!(tree::BnBTreeObj,node::BnBNode)
    ps = tree.options.log_levels
    BnBTree.check_print(ps,[:All,:FuncCall]) && println("update_incumbent")

    factor = tree.obj_fac
    if tree.incumbent == nothing || factor*node.best_bound > factor*tree.incumbent.objval
        objval = node.best_bound
        solution = copy(node.solution)
        status = :Optimal
        tree.incumbent = IncumbentSolution(objval,solution,status,tree.best_bound)
        if !tree.options.all_solutions 
            bound!(tree)
        end
        return true
    end
 
    return false
end


"""
    bound!(tree::BnBTreeObj)
"""
function bound!(tree::BnBTreeObj)
    function isbetter(n)
        return f*n.best_bound >= f*incumbent_val
    end
    incumbent_val = tree.incumbent.objval
    f = tree.obj_fac
    filter!(isbetter,tree.branch_nodes)
end

function print_table_header(fields, field_chars)
    ln = ""
    i = 1
    for f in fields
        padding = field_chars[i]-length(f)
        ln *= repeat(" ",trunc(Int, floor(padding/2)))
        ln *= f
        ln *= repeat(" ",trunc(Int, ceil(padding/2)))
        i += 1
    end
    println(ln)
    println(repeat("=", sum(field_chars)))
end

function is_table_diff(fields,last_arr,new_arr)
    if length(last_arr) != length(new_arr)
        return true
    end    

    time_idx = findfirst(fields .== "Time")
    if time_idx != 0
        last_arr = vcat(last_arr[1:time_idx-1],last_arr[time_idx+1:end])
        new_arr = vcat(new_arr[1:time_idx-1],new_arr[time_idx+1:end])
    end
    for i=1:length(last_arr)
        last_arr[i] != new_arr[i] && return true
    end
    return false 
end

function print_table(tree,node,step_obj,start_time,fields,field_chars,gain_gap;last_arr=[])
    arr = []
    nrestarts = step_obj.nrestarts
    i = 1
    ln = ""
    for f in fields
        val = ""
        if f == "Incumbent"
            val = tree.incumbent != nothing ? string(round(tree.incumbent.objval,2)) : "-"
        elseif f == "Best Bound"
            val = string(round(tree.best_bound,2))
        elseif f == "Gap"
            if tree.incumbent != nothing
                b = tree.best_bound
                f = tree.incumbent.objval
                val = string(round(abs(b-f)/abs(f)*100,1))*"%"
            else
                val = "-"
            end
        elseif f == "Time"
            val = string(round(time()-start_time,1))
        elseif f == "#Restarts"
            if nrestarts == -1
                val = "-"
            else
                val = string(nrestarts)
            end
        elseif f == "CLevel"
            val = string(node.level+1)
        elseif f == "GainGap"
            if gain_gap == Inf
                val = "∞"
            elseif gain_gap == -1.0
                val = "-"
            else
                val = string(round(gain_gap))*"%"
            end
            if length(val) > field_chars[i]
                val = ">>"
            end
        end
        padding = field_chars[i]-length(val)
        ln *= repeat(" ",trunc(Int, floor(padding/2)))
        ln *= val
        ln *= repeat(" ",trunc(Int, ceil(padding/2)))
        push!(arr,val)
        i += 1
    end
    
    BnBTree.is_table_diff(fields, last_arr,arr) && println(ln)
    return arr
end


function get_table_config(opts)
    if opts.branch_strategy == :StrongPseudoCost && opts.strong_restart
        fields = ["CLevel","Incumbent","Best Bound","Gap","Time","#Restarts"]
        field_chars = [8,28,28,7,8,10]
    else
        fields = ["CLevel","Incumbent","Best Bound","Gap","Time"]
        field_chars = [8,28,28,7,8]
    end
    
    if opts.branch_strategy == :StrongPseudoCost || opts.branch_strategy == :PseudoCost
        push!(fields, "GainGap")
        push!(field_chars, 10)
    end
    return fields, field_chars
end

function break_new_incumbent_limits(tree)
    if !isnan(tree.options.best_obj_stop)
        inc_val = tree.incumbent.objval
        bos = tree.options.best_obj_stop
        sense = tree.m.obj_sense
        if (sense == :Min && inc_val <= bos) || (sense == :Max && inc_val >= bos) 
            incu = tree.incumbent
            tree.incumbent = IncumbentSolution(incu.objval,incu.solution,:UserLimit,tree.best_bound)
            return true
        end
    end

    if tree.options.mip_gap != 0
        b = tree.best_bound
        f = tree.incumbent.objval
        gap_perc = abs(b-f)/abs(f)*100
        if gap_perc <= tree.options.mip_gap
            incu = tree.incumbent
            tree.incumbent = IncumbentSolution(incu.objval,incu.solution,:UserLimit,tree.best_bound)
            return true
        end
    end
    return false
end

function break_time_limit(tree)
    if !isnan(tree.options.time_limit) && time()-tree.start_time >= tree.options.time_limit
        if tree.incumbent == nothing
            tree.incumbent = IncumbentSolution(NaN,zeros(tree.m.num_var),:UserLimit,tree.best_bound)
            return true
        else
            tree.incumbent = IncumbentSolution(tree.incumbent.objval,tree.incumbent.solution,:UserLimit,tree.best_bound)
            return true
        end
    end
    return false
end

function add_obj_constr(tree)
    # add constr for objval
    if tree.options.incumbent_constr
        obj_expr = MathProgBase.obj_expr(tree.m.d)
        if tree.m.obj_sense == :Min
            obj_constr = Expr(:call, :<=, obj_expr, tree.incumbent.objval)
        else
            obj_constr = Expr(:call, :>=, obj_expr, tree.incumbent.objval)
        end
        MINLPBnB.expr_dereferencing!(obj_constr, tree.m.model)            
        # TODO: Change RHS instead of adding new (doesn't work for NL constraints atm)    
        JuMP.addNLconstraint(tree.m.model, obj_constr)
    end
end

function print_branch_nodes(nodes)
    for node in nodes
        println("------")
        println("BB: ", node.best_bound)
        println("LV: ", node.l_var)
        println("UV: ", node.u_var)
        println("------")
    end
end



"""
    solve(tree::BnBTreeObj,temp)

Solve the MIP part of a problem given by BnBTreeObj using branch and bound.
 - Identify the node to branch on
 - Get variable to branch on
 - Solve subproblems
"""
function solve(tree::BnBTreeObj)
    time_solve_leaf!s_get_idx = 0.0
    time_solve_leaf!s_branch = 0.0
    time_upd_gains = 0.0
    time_get_idx = 0.0
    time_branch = 0.0
    time_solve_leaf!s = 0.0
    time_bnb_solve_start = time()

    ps = tree.options.log_levels

    # check if already integral
    if BnBTree.are_type_correct(tree.m.solution,tree.m.var_type)
        tree.nsolutions = 1
        return tree.m
    end

    # Print table init
    if BnBTree.check_print(ps,[:Table]) 
        fields, field_chars = get_table_config(tree.options)
        print_table_header(fields,field_chars)
        last_table_arr = []
    end
    
    
    BnBTree.check_print(ps,[:All,:FuncCall]) && println("Solve Tree")
    
    # get variable where to split
    counter = 1    
    branch_strat = tree.options.branch_strategy

        
    first_incumbent = true
    gain_gap = 0
    
    all_nodes = []
    while length(tree.branch_nodes) > 0
        node = shift!(tree.branch_nodes)
        push!(all_nodes,node)

        step_obj = BnBTree.new_default_step_obj(node)
        
        BnBTree.upd_int_variable_idx!(step_obj,tree,counter)
        if step_obj.state == :Infeasible
            tree.incumbent = IncumbentSolution(NaN,zeros(tree.m.num_var),:Infeasible, NaN)
            break
        end

        BnBTree.branch!(tree,step_obj,counter)
        step_obj.state == :Break && break
        length(tree.branch_nodes) == 0 && break

        tree.best_bound = tree.obj_fac*maximum([tree.obj_fac*n.best_bound for n in tree.branch_nodes])

        # check if best
        if !tree.options.all_solutions && tree.incumbent != nothing && tree.incumbent.objval == tree.best_bound
            break
        end
        
        # maybe break on solution_limit (can be higher if two solutions found in last step)
        if tree.options.solution_limit > 0 && tree.nsolutions >= tree.options.solution_limit
            incu = tree.incumbent
            tree.incumbent = IncumbentSolution(incu.objval,incu.solution,:UserLimit,tree.best_bound)
            break
        end
    
        if BnBTree.check_print(ps,[:Table]) 
            if tree.options.branch_strategy != :StrongPseudoCost || counter > tree.options.strong_branching_nsteps
                step_obj.nrestarts = -1 # will be displayed as -
            end
            if counter <= tree.options.strong_branching_nsteps
                gain_gap = -1.0 # will be displayed as -
            end
            last_table_arr = print_table(tree,node,step_obj,time_bnb_solve_start,fields,field_chars,gain_gap;last_arr=last_table_arr)
        end

        if BnBTree.break_time_limit(tree)
            break
        end
        counter += 1
        time_solve_leaf!s_get_idx += step_obj.leaf_idx_time
        time_solve_leaf!s_branch += step_obj.leaf_branch_time
        time_branch += step_obj.branch_time
        time_get_idx += step_obj.idx_time
        time_upd_gains += step_obj.upd_gains_time
    end

    if tree.incumbent == nothing
        # infeasible
        tree.incumbent = IncumbentSolution(NaN, zeros(tree.m.num_var), :Infeasible, tree.best_bound)
    end

    # update best bound in incumbent
    tree.incumbent.best_bound = tree.best_bound

    if !isnan(tree.options.best_obj_stop)
        inc_val = tree.incumbent.objval
        bos = tree.options.best_obj_stop
        sense = tree.m.obj_sense
        if (sense == :Min && inc_val > bos) || (sense == :Max && inc_val < bos)
            warn("best_obj_gap couldn't be reached.")
        end
    end
    
    # println("ALL NODES")
    # print_branch_nodes(all_nodes)

    println("Incumbent status: ", tree.incumbent.status)

    time_bnb_solve = time()-time_bnb_solve_start
    println("#branches: ", counter)
    println("BnB time: ", round(time_bnb_solve,2))
    println("Solve leaf time get idx: ", round(time_solve_leaf!s_get_idx,2))
    println("Solve leaf time branch: ", round(time_solve_leaf!s_branch,2))
    println("Branch time: ", round(time_branch,2))
    println("Get idx time: ", round(time_get_idx,2))
    println("Upd gains time: ", round(time_upd_gains,2))
    return tree.incumbent
end

end
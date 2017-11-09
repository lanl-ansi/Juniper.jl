include("table_log.jl")

rtol = 1e-6
atol = 1e-6

type BnBNode
    idx                 :: Int64
    level               :: Int64
    l_var               :: Vector{Float64}
    u_var               :: Vector{Float64}
    solution            :: Vector{Float64}
    var_idx             :: Int64
    state               :: Symbol
    relaxation_state    :: Symbol
    best_bound          :: Union{Void,Float64}
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
    branch_nodes:: Vector{BnBNode}
    best_bound  :: Float64
end

# the object holds information for the current step
type StepObj
    node        :: BnBNode # current branch node
    var_idx     :: Int64   # variable to branch on
    state       :: Symbol  # if infeasible => break (might be set by strong branching)
    nrestarts   :: Int64 
    gain_gap    :: Float64
    gains       :: Vector{Float64}
    strong_int_vars :: Vector{Int64}
    idx_time    :: Float64
    leaf_idx_time :: Float64
    upd_gains_time :: Float64
    leaf_branch_time :: Float64
    branch_time     :: Float64
    integral       :: Vector{BnBNode}
    branch       :: Vector{BnBNode}
    l_nd        :: Union{Void,BnBNode}
    r_nd        :: Union{Void,BnBNode}
    counter     :: Int64
end

type TimeObj
    solve_leaves_get_idx :: Float64
    solve_leaves_branch :: Float64
    branch :: Float64
    get_idx :: Float64
    upd_gains :: Float64
end

include("bb_inits_and_defaults.jl")
include("bb_strategies.jl")
include("bb_user_limits.jl")
include("bb_type_correct.jl")
include("bb_integral_or_branch.jl")
include("bb_gains.jl")

function check_print(vec::Vector{Symbol}, ps::Vector{Symbol})
    for v in vec
        if v in ps
            return true
        end
    end
    return false
end

"""
    upd_int_variable_idx!(m,step_obj,opts,node,counter=1)    

Get the index of a variable to branch on.
"""
function upd_int_variable_idx!(m,step_obj,opts,int2var_idx,gain,gain_c,counter::Int64=1)  
    start = time()
    node = step_obj.node
    idx = 0
    strong_restarts = 0
    branch_strat = opts.branch_strategy
    status = :Normal
    if branch_strat == :MostInfeasible
        idx = branch_mostinfeasible(m,node)
    elseif branch_strat == :PseudoCost || branch_strat == :StrongPseudoCost
        if counter == 1 && branch_strat == :PseudoCost
            idx = branch_mostinfeasible(m,node)
        elseif counter <= opts.strong_branching_nsteps && branch_strat == :StrongPseudoCost
            status, idx, strong_restarts = branch_strong!(m,opts,int2var_idx,step_obj,counter)
        else
            idx = branch_pseudo(m,node,int2var_idx,gain,gain_c)
        end
    end
    step_obj.state = status
    step_obj.var_idx = idx
    step_obj.nrestarts = strong_restarts
    step_obj.idx_time = time()-start
    return
end

"""
    solve_leaf!(tree,leaf)

Solve a leaf by relaxation leaf is just a node.
Set the state and best_bound property
Update incumbent if new and add node to branch list if :Branch
Return state
"""
function solve_leaf!(m,step_obj,leaf,temp)
     # set bounds
    for i=1:m.num_var
        JuMP.setlowerbound(m.x[i], leaf.l_var[i])    
        JuMP.setupperbound(m.x[i], leaf.u_var[i])
    end

    status = JuMP.solve(m.model)
    objval = getobjectivevalue(m.model)
    leaf.solution = getvalue(m.x)
    status = status
    leaf.relaxation_state = status
    if status == :Error
        leaf.state = :Error
    elseif status == :Optimal
        leaf.best_bound = objval
        push_integral_or_branch!(m,step_obj,leaf,temp)
    else
        leaf.state = :Infeasible
    end
    return leaf.state
end

"""
    branch!(tree,step_obj,counter)

Branch a node by using x[idx] <= floor(x[idx]) and x[idx] >= ceil(x[idx])
Solve both nodes and set current node state to done.
"""
function branch!(m,opts,step_obj,counter;temp=false)
    ps = opts.log_levels
    node = step_obj.node
    vidx = step_obj.var_idx

    start = time()
    # it might be already branched on
    if node.state != :Branch
        for leaf in [step_obj.l_nd,step_obj.r_nd]
            if leaf.state == :Branch || leaf.state == :Integral
                push_integral_or_branch!(m,step_obj,leaf,false)
            end
        end
        return step_obj.l_nd,step_obj.r_nd
    end
    
    l_nd = new_default_node(node.idx*2,node.level+1,node.l_var,node.u_var,node.solution)
    r_nd = new_default_node(node.idx*2+1,node.level+1,node.l_var,node.u_var,node.solution)

    l_nd.u_var[vidx] = floor(node.solution[vidx])
    r_nd.l_var[vidx] = ceil(node.solution[vidx])

    # save that this node branches on this particular variable
    node.var_idx = vidx

    check_print(ps,[:All,:FuncCall]) && println("branch")
    
    if !temp
        node.state = :Done
        step_obj.l_nd = l_nd
        step_obj.r_nd = r_nd
    end

    start_leaf = time()
    l_state = solve_leaf!(m,step_obj,l_nd,temp)
    r_state = solve_leaf!(m,step_obj,r_nd,temp)
    leaf_time = time() - start_leaf

    if temp
        step_obj.leaf_idx_time += leaf_time
    else
        step_obj.leaf_branch_time += leaf_time
    end

    if l_state == :Break || r_state == :Break
        step_obj.state = :Break
    end

    branch_strat = opts.branch_strategy

    if check_print(ps,[:All])
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
    update_incumbent!(tree::BnBTreeObj,node::BnBNode)

Get's called if new integral solution was found. 
Check whether it's a new incumbent and update if necessary
"""
function update_incumbent!(tree::BnBTreeObj,node::BnBNode)
    ps = tree.options.log_levels
    check_print(ps,[:All,:FuncCall]) && println("update_incumbent")

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
        return f*n.best_bound > f*incumbent_val
    end
    incumbent_val = tree.incumbent.objval
    f = tree.obj_fac
    filter!(isbetter,tree.branch_nodes)
end

"""
    add_incumbent_constr(tree)

Add a constraint >=/<= incumbent 
"""
function add_incumbent_constr(tree)
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
        tree.m.ncuts += 1
    end
end

"""
    add_obj_epsilon_constr(tree)

Add a constraint obj ≦ (1+ϵ)*LB or obj ≧ (1-ϵ)*UB
"""
function add_obj_epsilon_constr(tree)
    # add constr for objval
    if tree.options.obj_epsilon > 0
        ϵ = tree.options.obj_epsilon
        obj_expr = MathProgBase.obj_expr(tree.m.d)
        if tree.m.obj_sense == :Min
            obj_constr = Expr(:call, :<=, obj_expr, (1+ϵ)*tree.m.objval)
        else
            obj_constr = Expr(:call, :>=, obj_expr, (1-ϵ)*tree.m.objval)
        end
        MINLPBnB.expr_dereferencing!(obj_constr, tree.m.model)            
        JuMP.addNLconstraint(tree.m.model, obj_constr)
        tree.m.ncuts += 1
    end
end

"""
    get_next_branch_node!(tree)

Get the next branch node (sorted by best bound)
Return true,step_obj if there is a branch node and
false, nothing otherwise
"""
function get_next_branch_node!(tree)
    if length(tree.branch_nodes) == 0
        return false,nothing
    end
    bvalue, nidx = findmax([tree.obj_fac*n.best_bound for n in tree.branch_nodes])

    trav_strat = tree.options.traverse_strategy
    if trav_strat == :DFS || (trav_strat == :DBFS && tree.incumbent == nothing)
        _value, nidx = findmax([n.level for n in tree.branch_nodes])
    end

    node = tree.branch_nodes[nidx]
    deleteat!(tree.branch_nodes,nidx)

    tree.best_bound = tree.obj_fac*bvalue
    bbreak = isbreak_mip_gap(tree)
    if bbreak 
        return false,nothing
    end
    return true,new_default_step_obj(tree.m,node)
end


"""
    one_branch_step!(m, opts, step_obj,int2var_idx,gain,gain_c, counter)

Get a branch variable using the specified strategy and branch on the node in step_obj 
using that variable. Return the new updated step_obj
"""
function one_branch_step!(m, opts, step_obj,int2var_idx,gain,gain_c, counter)
    node = step_obj.node
    step_obj.counter = counter

# get branch variable    
    upd_int_variable_idx!(m,step_obj,opts,int2var_idx,gain,gain_c,counter)
    if step_obj.state != :GlobalInfeasible && step_obj.state != :LocalInfeasible
        # branch
        branch!(m,opts,step_obj,counter)
    end
    return step_obj
end

"""
    upd_time_obj!(time_obj, step_obj)

Add step_obj times to time_obj
"""
function upd_time_obj!(time_obj, step_obj)
    time_obj.solve_leaves_get_idx += step_obj.leaf_idx_time
    time_obj.solve_leaves_branch += step_obj.leaf_branch_time
    time_obj.branch += step_obj.branch_time
    time_obj.get_idx += step_obj.idx_time
    time_obj.upd_gains += step_obj.upd_gains_time
end

function init_time_obj()
    return TimeObj(0.0,0.0,0.0,0.0,0.0)
end

"""
    upd_tree_obj!(tree,step_obj,time_obj)

Update the tree obj like new incumbent or new branch nodes using the step_obj
Return false if it's the end of the algorithm (checking different break rules)
"""
function upd_tree_obj!(tree,step_obj,time_obj)
    node = step_obj.node
    still_running = true

    if step_obj.node.level+1 > tree.m.nlevels
        tree.m.nlevels = step_obj.node.level+1
    end

    if step_obj.state == :GlobalInfeasible
        tree.incumbent = IncumbentSolution(NaN,zeros(tree.m.num_var),:Infeasible, NaN)
        still_running = false 
    end

    if step_obj.state == :Break 
        still_running = false
    end
    
    if still_running
        upd_gains_step!(tree,step_obj)
    end

    bbreak = upd_integral_branch!(tree,step_obj)
    if bbreak 
        still_running = false 
    end

    upd_time_obj!(time_obj,step_obj)
    if still_running
        return false
    else
        return true
    end
end

"""
    solve_sequential(tree,
        last_table_arr,
        time_bnb_solve_start,
        fields,
        field_chars,
        time_obj)
    
Run branch and bound on a single processor
"""
function solve_sequential(tree,
    last_table_arr,
    time_bnb_solve_start,
    fields,
    field_chars,
    time_obj)

    m = tree.m
    opts = tree.options
    int2var_idx = tree.int2var_idx
    gain = tree.obj_gain
    gain_c = tree.obj_gain_c
    counter = 1
    ps = tree.options.log_levels
    while true
        exists,step_obj = get_next_branch_node!(tree)
        !exists && break
        isbreak_after_step!(tree) && break
        step_obj = one_branch_step!(m, opts, step_obj,int2var_idx,gain,gain_c, counter)
        m.nnodes += 2 # two nodes explored per branch
        node = step_obj.node

        bbreak = upd_tree_obj!(tree,step_obj,time_obj)
        
        if check_print(ps,[:Table]) 
            last_table_arr = print_table(1,tree,node,step_obj,time_bnb_solve_start,fields,field_chars,counter;last_arr=last_table_arr)
        end

        if bbreak 
            break
        end
        counter += 1
    end
    return counter
end

"""
    pmap(f, tree, counter, last_table_arr, time_bnb_solve_start,
        fields, field_chars, time_obj)

Run the solving steps on several processors
"""
function pmap(f, tree, counter, last_table_arr, time_bnb_solve_start,
    fields, field_chars, time_obj)
    np = nprocs()  # determine the number of processes available
    if np < tree.options.processors
        warn("Julia was started with less processors then you define in your options")
    end
    if tree.options.processors < np
        np = tree.options.processors
    end

    # function to produce the next work item from the queue.
    # in this case it's just an index.
    ps = tree.options.log_levels
    still_running = true
    run_counter = 0
    counter = 0

    for p=1:np
        remotecall_fetch(srand, p, 1)
    end

    branch_strat = tree.options.branch_strategy
    opts = tree.options
    @sync begin
        for p=1:np
            if p != myid() || np == 1
                @async begin
                    while true
                        exists,step_obj = get_next_branch_node!(tree)

                        while !exists && still_running
                            sleep(0.1)
                            exists,step_obj = get_next_branch_node!(tree)
                            exists && break
                        end
                        if !still_running
                            break
                        end
                        
                        if isbreak_after_step!(tree) 
                            println("is break after step")
                            still_running = false 
                            break
                        end
                        
                        counter += 1
                        run_counter += 1
                        step_obj = remotecall_fetch(f, p, tree.m, tree.options, step_obj, tree.int2var_idx,tree.obj_gain,tree.obj_gain_c, counter)
                        tree.m.nnodes += 2 # two nodes explored per branch
                        run_counter -= 1
                        !still_running && break
                    
                        bbreak = upd_tree_obj!(tree,step_obj,time_obj)

                        if run_counter == 0 && length(tree.branch_nodes) == 0
                            still_running = false 
                        end

                        if bbreak
                            still_running = false
                        end

                        if check_print(ps,[:Table]) 
                            last_table_arr = print_table(p,tree,step_obj.node,step_obj,time_bnb_solve_start,fields,field_chars,counter;last_arr=last_table_arr)
                        end

                        if !still_running
                            break
                        end
                    end
                end
            end
        end
    end
    return counter
end

"""
    solvemip(tree::BnBTreeObj)

Solve the MIP part of a problem given by BnBTreeObj using branch and bound.
 - Identify the node to branch on
 - Get variable to branch on
 - Solve subproblems
"""
function solvemip(tree::BnBTreeObj)
    time_obj = init_time_obj()
    time_bnb_solve_start = time()

    ps = tree.options.log_levels

    # check if already integral
    if are_type_correct(tree.m.solution,tree.m.var_type)
        tree.nsolutions = 1
        return tree.m
    end

    last_table_arr = []
    fields = []
    field_chars = []
    # Print table init
    if check_print(ps,[:Table]) 
        fields, field_chars = get_table_config(tree.options)
        print_table_header(fields,field_chars)        
    end
    
    
    check_print(ps,[:All,:FuncCall]) && println("Solve Tree")
    
    counter = 1    
    branch_strat = tree.options.branch_strategy

    add_obj_epsilon_constr(tree)

    # use pmap if more then one processor
    if tree.options.processors > 1
        counter = pmap(MINLPBnB.one_branch_step!,tree,
            counter,
            last_table_arr,
            time_bnb_solve_start,
            fields,
            field_chars,
            time_obj)
    else 
        counter = solve_sequential(tree,
            last_table_arr,
            time_bnb_solve_start,
            fields,
            field_chars,
            time_obj)
    end
    counter -= 1 # to have the correct number of branches

    if tree.incumbent == nothing
        # infeasible
        tree.incumbent = IncumbentSolution(NaN, zeros(tree.m.num_var), :Infeasible, tree.best_bound)
    end

    # update best bound in incumbent
    tree.incumbent.best_bound = tree.best_bound

    if tree.options.obj_epsilon != 0 && tree.incumbent.status == :Infeasible
        warn("Maybe only infeasible because of obj_epsilon.")
    end

    if !isnan(tree.options.best_obj_stop)
        inc_val = tree.incumbent.objval
        bos = tree.options.best_obj_stop
        sense = tree.m.obj_sense
        if (sense == :Min && inc_val > bos) || (sense == :Max && inc_val < bos)
            warn("best_obj_gap couldn't be reached.")
        end
    end
    
    tree.m.nbranches = counter

    time_bnb_solve = time()-time_bnb_solve_start
    (:Table in tree.options.log_levels) && println("")
    (:Info in tree.options.log_levels) && println("#branches: ", counter)
    if :Timing in tree.options.log_levels
        println("BnB time: ", round(time_bnb_solve,2))
        println("% leaf time: ", round((time_obj.solve_leaves_get_idx+time_obj.solve_leaves_branch)/time_bnb_solve*100,1))
        println("Solve leaf time get idx: ", round(time_obj.solve_leaves_get_idx,2))
        println("Solve leaf time branch: ", round(time_obj.solve_leaves_branch,2))
        println("Branch time: ", round(time_obj.branch,2))
        println("Get idx time: ", round(time_obj.get_idx,2))
        println("Upd gains time: ", round(time_obj.upd_gains,2))
    end
    return tree.incumbent
end
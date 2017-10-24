type MINLPBnBModel <: MathProgBase.AbstractNonlinearModel
    nl_solver       :: MathProgBase.AbstractMathProgSolver
    print_syms      :: Vector{Symbol}
    model           :: Union{Void,JuMP.Model}
        
    status          :: Symbol
    objval          :: Float64

    x               
    num_constr      :: Int64
    num_nl_constr   :: Int64
    num_l_constr    :: Int64
    num_var         :: Int64
    obj_expr        
    l_var           :: Vector{Float64}
    u_var           :: Vector{Float64}
    l_constr        :: Vector{Float64}
    u_constr        :: Vector{Float64}

    var_type        :: Vector{Symbol}
    constr_type     :: Vector{Symbol}
    isconstrlinear  :: Vector{Bool}
    obj_sense       :: Symbol
    d               :: MathProgBase.AbstractNLPEvaluator
    num_int_bin_var :: Int64

    solution        :: Vector{Float64}

    soltime         :: Float64

    MINLPBnBModel() = new()
end

"""
    MathProgBase.NonlinearModel(s::MINLPBnBSolverObj)

Generate NonLinearModel and specify nl solver
"""
function MathProgBase.NonlinearModel(s::MINLPBnBSolverObj)
    println("model.jl => MathProgBase.NonlinearModel")
    return MINLPBnBNonlinearModel(s.nl_solver,s.print_syms)
end

"""
    MINLPBnBNonlinearModel(lqps::MathProgBase.AbstractMathProgSolver)

Initialize the NonLinearModel with the solver, set status, objval and solution
"""
function MINLPBnBNonlinearModel(lqps::MathProgBase.AbstractMathProgSolver,print_syms)
    println("model.jl => MINLPBnBNonlinearModel")
    m = MINLPBnBModel() # don't initialise everything yet

    m.nl_solver = lqps
    m.print_syms = print_syms
    m.status = :None
    m.objval = NaN
    m.solution = Float64[]

    return m
end

"""
    MathProgBase.loadproblem!(m,num_var,num_constr,l_var,u_var,l_constr,u_constr,sense,d)

Initialize other fields MINLPBnBModel after all variables, constraints and the objective is set
"""
function MathProgBase.loadproblem!(
    m::MINLPBnBModel,
    num_var::Int, num_constr::Int,
    l_var::Vector{Float64}, u_var::Vector{Float64},
    l_constr::Vector{Float64}, u_constr::Vector{Float64},
    sense::Symbol, d::MathProgBase.AbstractNLPEvaluator)

    # initialise other fields
    m.num_var = num_var
    m.num_constr = num_constr
    m.l_var    = l_var
    m.u_var    = u_var
    m.l_constr = l_constr
    m.u_constr = u_constr
    m.d = d
    m.obj_sense = sense
    m.solution = fill(NaN, m.num_var)
    m.var_type = fill(:Cont,num_var)

    println("loadproblem!")
    MathProgBase.initialize(m.d, [:Grad,:Jac,:Hess,:ExprGraph])
end

#=
    Used from https://github.com/lanl-ansi/POD.jl
=# 
function expr_dereferencing(expr, m)
    for i in 2:length(expr.args)
        if isa(expr.args[i], Float64)
            k = 0
        elseif expr.args[i].head == :ref
            @assert isa(expr.args[i].args[2], Int)
            expr.args[i] = Variable(m, expr.args[i].args[2])
        elseif expr.args[i].head == :call
            expr_dereferencing(expr.args[i], m)
        else
            error("expr_dereferencing :: Unexpected term in expression tree.")
        end
    end
end

"""
    divide_nl_l_constr(m::MINLPBnBModel)

Get # of linear and non linear constraints and save for each index if linear or non linear    
"""
function divide_nl_l_constr(m::MINLPBnBModel)
    isconstrlinear = Array{Bool}(m.num_constr)
    m.num_l_constr = 0
    for i = 1:m.num_constr
        isconstrlinear[i] = MathProgBase.isconstrlinear(m.d, i)
        if isconstrlinear[i]
            m.num_l_constr += 1
        end
    end
    m.num_nl_constr = m.num_constr - m.num_l_constr  
    m.isconstrlinear = isconstrlinear
end

function replace_solution!(old_m::MINLPBnBModel, new_m::MINLPBnBModel)
    old_m.solution = new_m.solution
    old_m.objval = new_m.objval
    old_m.model = new_m.model    
end

"""
    MathProgBase.optimize!(m::MINLPBnBModel)

Optimize by creating a model based on the variables saved in MINLPBnBModel.
"""
function MathProgBase.optimize!(m::MINLPBnBModel)
    println("optimize!")
    println("Types")
    println(m.var_type)

    m.model = Model(solver=m.nl_solver) 
    lb = [m.l_var; -1e6]
    ub = [m.u_var; 1e6]
    # all continuous 
    @variable(m.model, lb[i] <= x[i=1:m.num_var] <= ub[i])

    # define the objective function
    obj_expr = MathProgBase.obj_expr(m.d)
    expr_dereferencing(obj_expr, m.model)
    JuMP.setNLobjective(m.model, m.obj_sense,  obj_expr)

    divide_nl_l_constr(m)

    # add all constraints
    for i=1:m.num_constr
        constr_expr = MathProgBase.constr_expr(m.d,i)
        expr_dereferencing(constr_expr, m.model)
        JuMP.addNLconstraint(m.model, constr_expr)
    end

    m.x = x
    start = time()
    status = solve(m.model)
    m.soltime = time()-start
    println("Time for relaxation: ", m.soltime)

    m.objval   = getobjectivevalue(m.model)
    m.solution = getvalue(x)

    println("Solution: ", m.solution)

    m.status = status

    bnbtree = BnBTree.init(m)
    bnbtree_m = BnBTree.solve(bnbtree)

    replace_solution!(m, bnbtree_m)
    m.soltime = time()-start
    return m.status

end

MathProgBase.setwarmstart!(m::MINLPBnBModel, x) = fill(0.0, length(x))

"""
    MathProgBase.setvartype!(m::MINLPBnBModel, v::Vector{Symbol}) 

Is called between loadproblem! and optimize! and has a vector v of types for each variable.
The number of int/bin variables is saved in num_int_bin_var
"""
function MathProgBase.setvartype!(m::MINLPBnBModel, v::Vector{Symbol}) 
    m.var_type = v
    c = count(i->(i==:Int || i==:Bin), v)
    m.num_int_bin_var = c
end

MathProgBase.status(m::MINLPBnBModel) = m.status
MathProgBase.getobjval(m::MINLPBnBModel) = getobjectivevalue(m.model)

# any auxiliary variables will need to be filtered from this at some point
MathProgBase.getsolution(m::MINLPBnBModel) = m.solution

MathProgBase.getsolvetime(m::MINLPBnBModel) = m.soltime
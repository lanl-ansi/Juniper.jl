type MINLPBnBModel <: MathProgBase.AbstractNonlinearModel
    nl_solver       :: MathProgBase.AbstractMathProgSolver
    model           :: Union{Void,JuMP.Model}
        
    status          :: Symbol
    objval          :: Float64

    num_constr      :: Int64
    num_nl_constr   :: Int64
    num_l_constr    :: Int64
    num_var         :: Int64
    l_var           :: Vector{Float64}
    u_var           :: Vector{Float64}
    l_constr        :: Vector{Float64}
    u_constr        :: Vector{Float64}
    isobjlinear     :: Bool

    A
    A_lb
    A_ub
    l               :: Vector{Float64}
    u               :: Vector{Float64}
    c               :: Vector{Float64}
    var_type        :: Vector{Symbol}
    constr_type     :: Vector{Symbol}
    isconstrlinear  :: Vector{Bool}
    obj_sense       :: Symbol
    mip_x
    d               :: MathProgBase.AbstractNLPEvaluator
    num_int_var     :: Int64

    solution        :: Vector{Float64}

    soltime         :: Float64

    MINLPBnBModel() = new()
end

function MathProgBase.NonlinearModel(s::MINLPBnBSolverObj)
    println("model.jl => MathProgBase.NonlinearModel")
    return MINLPBnBNonlinearModel(s.nl_solver)
end

function MINLPBnBNonlinearModel(lqps::MathProgBase.AbstractMathProgSolver)
    println("model.jl => MINLPBnBNonlinearModel")
    m = MINLPBnBModel() # don't initialise everything yet

    m.nl_solver = lqps
    m.status = :None
    m.objval = NaN
    m.solution = Float64[]

    return m
end

function gen_linear(m::MINLPBnBModel)
    # set up map of linear rows
    isconstrlinear = Array{Bool}(m.num_constr)
    numlinear = 0
    constraint_to_linear = fill(-1,m.num_constr)
    for i = 1:m.num_constr
        isconstrlinear[i] = MathProgBase.isconstrlinear(m.d, i)
        if isconstrlinear[i]
            numlinear += 1
            constraint_to_linear[i] = numlinear
        end
    end
    m.num_nl_constr = m.num_constr - numlinear

    println("#linear constr: ", numlinear)

    # extract sparse jacobian structure
    jac_I, jac_J = MathProgBase.jac_structure(m.d)

    # evaluate jacobian at x = 0
    c = zeros(m.num_var)
    x = m.solution
    jac_V = zeros(length(jac_I))
    MathProgBase.eval_jac_g(m.d, jac_V, x)
    MathProgBase.eval_grad_f(m.d, c, x)
    m.isobjlinear = MathProgBase.isobjlinear(m.d)
    if m.isobjlinear
        println("Objective function is linear")
        m.c = c
    else
        println("Objective function is nonlinear")
        m.c = zeros(m.num_var)
    end

    # Build up sparse matrix for linear constraints
    A_I = Int[]
    A_J = Int[]
    A_V = Float64[]

    for k in 1:length(jac_I)
        row = jac_I[k]
        if !isconstrlinear[row]
            continue
        end
        row = constraint_to_linear[row]
        push!(A_I,row); push!(A_J, jac_J[k]); push!(A_V, jac_V[k])
    end

    m.A = sparse(A_I, A_J, A_V, numlinear, m.num_var)

    # g(x) might have a constant, i.e., a'x + b
    # let's find b
    constraint_value = zeros(m.num_constr)
    MathProgBase.eval_g(m.d, constraint_value, x)

    b = constraint_value[isconstrlinear] - m.A * x

    # so linear constraints are of the form lb ≤ a'x + b ≤ ub

    # set up A_lb and A_ub vectors
    m.A_lb = m.l_constr[isconstrlinear] - b
    m.A_ub = m.u_constr[isconstrlinear] - b

    # Now we have linear parts
    m.isconstrlinear = isconstrlinear
end

function gen_linear_model(m::MINLPBnBModel, model)
    lb = [m.l_var; -1e6]
    ub = [m.u_var; 1e6]
    @variable(model, lb[i] <= x[i=1:m.num_var+1] <= ub[i])
    num_int_var = 0
    for i = 1:m.num_var
        JuMP.setcategory(x[i], m.var_type[i])
        if m.var_type[i] == :Int || m.var_type[i] == :Bin
            num_int_var += 1
        end
    end
    if num_int_var == 0
        error("No variables of type integer or binary; call the conic continuous solver directly for pure continuous problems")
    end
    JuMP.setcategory(x[m.num_var+1], :Cont)
    for i = 1:m.num_constr-m.num_nl_constr
        if m.A_lb[i] > -Inf && m.A_ub[i] < Inf
            if m.A_lb[i] == m.A_ub[i]
                @constraint(model, m.A[i:i,:]*x[1:m.num_var] .== m.A_lb[i])
            else
                @constraint(model, m.A[i:i,:]*x[1:m.num_var] .>= m.A_lb[i])
                @constraint(model, m.A[i:i,:]*x[1:m.num_var] .<= m.A_ub[i])
            end
        elseif m.A_lb[i] > -Inf
            @constraint(model, m.A[i:i,:]*x[1:m.num_var] .>= m.A_lb[i])
        else
            @constraint(model, m.A[i:i,:]*x[1:m.num_var] .<= m.A_ub[i])
        end
    end
    c_new = [m.obj_sense == :Max ? -m.c : m.c; m.isobjlinear ? 0.0 : 1.0]
    @objective(model, Min, dot(c_new, x))

    m.mip_x = x
    m.num_int_var = num_int_var
end

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
    println("typeof(m.d): ",typeof(m.d))

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

function MathProgBase.optimize!(m::MINLPBnBModel)
    println("optimize!")
    println("Types")
    println(m.var_type)

    # if we haven't gotten a starting point,
    # (e.g., if acting as MIQP solver) assume zero is okay
    if any(isnan,m.solution)
        m.solution = zeros(length(m.solution))
    end

    m.model = Model(solver=m.nl_solver) 
    gen_linear(m)
    gen_linear_model(m, m.model)

    for i=1:m.num_constr
        if !m.isconstrlinear[i]
            println("Constr #",i , " is non linear")
            constr_expr = MathProgBase.constr_expr(m.d,i)
            expr_dereferencing(constr_expr, m.model)
            JuMP.addNLconstraint(m.model, constr_expr)
            # @constraint(m.model, constr_expr)
        end
    end


    start = time()
    status = solve(m.model, relaxation=true)
    m.soltime = time()-start

    m.status = status
    return m.status

end

MathProgBase.setwarmstart!(m::MINLPBnBModel, x) = fill(0.0, length(x))

MathProgBase.setvartype!(m::MINLPBnBModel, v::Vector{Symbol}) = (m.var_type = v)

MathProgBase.status(m::MINLPBnBModel) = m.status
MathProgBase.getobjval(m::MINLPBnBModel) = getobjectivevalue(m.model)

# any auxiliary variables will need to be filtered from this at some point
MathProgBase.getsolution(m::MINLPBnBModel) = MathProgBase.getsolution(internalmodel(m.model))

MathProgBase.getsolvetime(m::MINLPBnBModel) = m.soltime

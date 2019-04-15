include("debug.jl")
include("fpump.jl")

"""
    MathProgBase.NonlinearModel(s::JuniperSolverObj)

Generate NonLinearModel and specify nl solver
"""
function MathProgBase.NonlinearModel(s::JuniperSolverObj)
    return JuniperNonlinearModel(s)
end

"""
    JuniperNonlinearModel(lqps::MathProgBase.AbstractMathProgSolver)

Initialize the NonLinearModel with the solver, set status, objval and solution
"""
function JuniperNonlinearModel(s::JuniperSolverObj)
    m = JuniperModel() # don't initialise everything yet

    m.nl_solver = s.nl_solver
    m.options = s.options
    m.status = :None
    m.objval = NaN
    m.best_bound = NaN
    m.solution = Float64[]
    m.nsolutions = 0
    m.solutions = []
    m.num_disc_var = 0
    m.nintvars = 0
    m.nbinvars = 0
    m.nnodes = 1 # is set to one for the root node
    m.ncuts = 0
    m.nbranches = 0
    m.nlevels = 1
    m.relaxation_time = 0.0
    if m.options.mip_solver != nothing
        m.mip_solver = m.options.mip_solver
    end

    return m
end

"""
    MathProgBase.loadproblem!(m,num_var,num_constr,l_var,u_var,l_constr,u_constr,sense,d)

Initialize other fields JuniperModel after all variables, constraints and the objective is set
"""
function MathProgBase.loadproblem!(
    m::JuniperModel,
    num_var::Int, num_constr::Int,
    l_var::Vector{Float64}, u_var::Vector{Float64},
    l_constr::Vector{Float64}, u_constr::Vector{Float64},
    sense::Symbol, d::MathProgBase.AbstractNLPEvaluator)

    VERSION > v"0.7.0-" ? Random.seed!(1) : srand(1)

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
    m.start_value = zeros(m.num_var)

    MathProgBase.initialize(m.d, [:ExprGraph,:Jac,:Grad])
end

function replace_solution!(m::JuniperModel, best_known)
    m.solution = best_known.solution
    m.objval = best_known.objval
    m.status = best_known.status
    m.best_bound = best_known.best_bound # is reasonable for gap or time limit
end

"""
    MathProgBase.optimize!(m::JuniperModel)

Optimize by creating a model based on the variables saved in JuniperModel.
"""
function MathProgBase.optimize!(m::JuniperModel)
    m.start_time = time()
    ps = m.options.log_levels
    m.debugDict = Dict{Any,Any}()
    if !m.options.fixed_gain_mu && m.obj_sense == :Max
        m.options.gain_mu = 1-m.options.gain_mu
    end
    (:All in ps || :AllOptions in ps) && print_options(m;all=true)
    (:Options in ps) && print_options(m;all=false)

    nw = nworkers()
    if nw < m.options.processors
        m.options.processors = nw
        @warn "Julia was started with less processors then you define in your options"
    end

    create_root_model!(m)
    relaxation_start_time = time()
    restarts = solve_root_model!(m)

    (:All in ps || :Info in ps) && println("Status of relaxation: ", m.status)

    m.soltime = time()-m.start_time
    m.relaxation_time = time()-relaxation_start_time
    m.options.debug && debug_fill_basic(m.debugDict,m,restarts)
    if m.status != :Optimal && m.status != :LocalOptimal
        if m.options.debug && m.options.debug_write
            write(m.options.debug_file_path, JSON.json(m.debugDict))
        end
        return m.status
    end

    (:All in ps || :Info in ps || :Timing in ps) && println("Time for relaxation: ", m.relaxation_time)
    m.objval   = getobjectivevalue(m.model)
    m.solution = getvalue(m.x)

    m.options.debug && debug_objective(m.debugDict,m)

    internal_model = internalmodel(m.model)
    if hasmethod(MathProgBase.freemodel!, Tuple{typeof(internal_model)})
        MathProgBase.freemodel!(internal_model)
    end

    (:All in ps || :Info in ps || :Timing in ps) && println("Relaxation Obj: ", m.objval)

    inc_sol, inc_obj = nothing, nothing
    if m.num_disc_var > 0
        if m.options.feasibility_pump
            inc_sol, inc_obj = fpump(m)
        end
        bnbtree = init(m.start_time, m; inc_sol = inc_sol, inc_obj = inc_obj)
        best_known = solvemip(bnbtree)

        replace_solution!(m, best_known)
        m.nsolutions = bnbtree.nsolutions
    else
        m.nsolutions = 1
        # TODO should be getobjbound but that is not working
        m.best_bound = getobjectivevalue(m.model)
    end
    m.soltime = time()-m.start_time

    (:All in ps || :Info in ps) && println("Obj: ",m.objval)

    if length(m.solutions) == 0
        push!(m.solutions, SolutionObj(m.solution, m.objval))
    end

    m.options.debug && debug_set_solution(m.debugDict,m)
    if m.options.debug && m.options.debug_write
        write(m.options.debug_file_path, JSON.json(m.debugDict))
    end
    return m.status
end

function create_root_model!(m::JuniperModel)
    ps = m.options.log_levels

    m.model = Model(solver=m.nl_solver)
    lb = m.l_var
    ub = m.u_var
    # all continuous we solve relaxation first
    @variable(m.model, lb[i] <= x[i=1:m.num_var] <= ub[i])
    setvalue(x[1:m.num_var], m.start_value)

    # define the objective function
    obj_expr = MathProgBase.obj_expr(m.d)
    expr_dereferencing!(obj_expr, m.model)
    JuMP.setNLobjective(m.model, m.obj_sense, obj_expr)

    divide_nl_l_constr(m)
    if m.num_l_constr > 0
        m.affs = construct_affine_vector(m)
    end
    (:All in ps || :Info in ps) && print_info(m)

    # add all constraints
    for i=1:m.num_constr
        if m.isconstrlinear[i]
            constr = m.affs[i]
            if constr.sense == :(>=)
                @constraint(m.model, sum(x[constr.var_idx[j]]*constr.coeff[j] for j=1:length(constr.var_idx)) >= constr.rhs)
            elseif constr.sense == :(<=)
                @constraint(m.model, sum(x[constr.var_idx[j]]*constr.coeff[j] for j=1:length(constr.var_idx)) <= constr.rhs)
            else # ==
                @constraint(m.model, sum(x[constr.var_idx[j]]*constr.coeff[j] for j=1:length(constr.var_idx)) == constr.rhs)
            end
        else
            constr_expr = MathProgBase.constr_expr(m.d,i)
            expr_dereferencing!(constr_expr, m.model)
            JuMP.addNLconstraint(m.model, constr_expr)
        end
    end

    m.x = x
end

function solve_root_model!(m::JuniperModel)
    m.status = solve(m.model; suppress_warnings=true)
    restarts = 0
    max_restarts = m.options.num_resolve_root_relaxation
    m.options.debug && debug_init(m.debugDict)
    while m.status != :Optimal && m.status != :LocalOptimal &&
        restarts < max_restarts && time()-m.start_time < m.options.time_limit

        internal_model = internalmodel(m.model)
        if hasmethod(MathProgBase.freemodel!, Tuple{typeof(internal_model)})
            MathProgBase.freemodel!(internal_model)
        end
        restart_values = generate_random_restart(m)
        m.options.debug && debug_restart_values(m.debugDict,restart_values)
        for i=1:m.num_var
            setvalue(m.x[i], restart_values[i])
        end
        m.status = solve(m.model)
        restarts += 1
    end
    return restarts
end

MathProgBase.setwarmstart!(m::JuniperModel, x) = (m.start_value = x)

"""
    MathProgBase.setvartype!(m::JuniperModel, v::Vector{Symbol})

Is called between loadproblem! and optimize! and has a vector v of types for each variable.
The number of int/bin variables is saved in num_disc_var
"""
function MathProgBase.setvartype!(m::JuniperModel, v::Vector{Symbol})
    m.var_type = v
    m.nintvars = count(i->(i==:Int), v)
    m.nbinvars = count(i->(i==:Bin), v)
    m.num_disc_var =  m.nintvars + m.nbinvars
    for (i,s) in enumerate(v)
        if s==:Bin
            m.l_var[i] = 0
            m.u_var[i] = 1
        end
    end
    m.disc2var_idx = zeros(m.num_disc_var)
    m.var2disc_idx = zeros(m.num_var)
    int_i = 1
    for i=1:m.num_var
        if m.var_type[i] != :Cont
            m.disc2var_idx[int_i] = i
            m.var2disc_idx[i] = int_i
            int_i += 1
        end
    end
end

MathProgBase.status(m::JuniperModel) = m.status

# any auxiliary variables will need to be filtered from this at some point
MathProgBase.getsolution(m::JuniperModel) = m.solution

MathProgBase.getsolvetime(m::JuniperModel) = m.soltime

MathProgBase.getobjval(m::JuniperModel) = m.objval

MathProgBase.getobjbound(m::JuniperModel) = m.best_bound

function MathProgBase.getobjgap(m::JuniperModel)
    b = m.best_bound
    if isnan(m.objval)
        return NaN
    else
        f = m.objval
        return abs(b-f)/abs(f)
    end
end

getnsolutions(m::JuniperModel) = m.nsolutions
getsolutions(m::JuniperModel) = m.solutions
getnbranches(m::JuniperModel) = m.nbranches
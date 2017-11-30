type Aff
    sense     :: Symbol
    var_idx   :: Vector{Int64}
    coeff     :: Vector{Float64}
    rhs       :: Float64

    Aff() = new()
end

type TabuList
    sols      :: Vector{Vector{Float64}}
    length    :: Int64
    pointer   :: Int64

    TabuList() = new()
end


function generate_mip(m, nlp_sol, tabu_list)
    mip_model = Model(solver=m.mip_solver)
    lb = [m.l_var; -1e6]
    ub = [m.u_var; 1e6]
    @variable(mip_model, lb[i] <= mx[i=1:m.num_var] <= ub[i])
    for i=1:m.num_var
        if m.var_type[i] == :Int
            setcategory(mx[i], :Int)
        elseif m.var_type[i] == :Bin
            setcategory(mx[i], :Bin)
        end
    end

    g = zeros(m.num_constr)
    MathProgBase.eval_g(m.d,g,nlp_sol)

    js = MathProgBase.jac_structure(m.d)

    jg = zeros(length(js[1]))
    MathProgBase.eval_jac_g(m.d, jg, nlp_sol)

    # Construct the data structure for our affine constraints
    aff = Vector{Aff}(m.num_constr)
    for i=1:m.num_constr
        aff[i] = Aff()
        aff[i].var_idx = []
        aff[i].coeff = []
        constr_expr = MathProgBase.constr_expr(m.d,i)
        aff[i].rhs = constr_expr.args[3]
        aff[i].sense = constr_expr.args[1]
    end

    # if linear constraint the derivative are the coeffs
    idx = 1
    for row in js[1]
        if m.isconstrlinear[row]
            col = js[2][idx]
            push!(aff[row].var_idx, col)
            push!(aff[row].coeff, jg[idx])
        end
        idx += 1
    end
    
    counter = 1
    linear = false
    constr_removed = 0
    for c in aff
        if m.isconstrlinear[counter]
            if c.sense == :(>=)
                @constraint(mip_model, sum(mx[c.var_idx[i]]*c.coeff[i] for i=1:length(c.var_idx)) >= c.rhs)
            elseif c.sense == :(<=)
                @constraint(mip_model, sum(mx[c.var_idx[i]]*c.coeff[i] for i=1:length(c.var_idx)) <= c.rhs)
            elseif c.sense == :(==) && linear
                @constraint(mip_model, sum(mx[c.var_idx[i]]*c.coeff[i] for i=1:length(c.var_idx)) == c.rhs)
            end
            counter += 1
        end
    end

    @variable(mip_model, mx_p[i=1:m.num_int_bin_var] >= 0)
    @variable(mip_model, mx_m[i=1:m.num_int_bin_var] >= 0)
    for i=1:m.num_int_bin_var
        vi = m.int2var_idx[i]
        @constraint(mip_model, mx[vi] == nlp_sol[vi]+mx_p[i]-mx_m[i])
    end
    
    #=
        For all discrete i
        x_{i} => mx[i]
        |x|   => m.num_int_bin_var
        x'_{ik} => solution for x[i] in tabu list (pointer -> k)
        j => the discrete index of i
        M(1-B_{jk}) ≧ x_{i}-x'_{ik}
        -M(1-B_{jk}) ≦ x_{i}-x'_{ik} 
        ∑ B ≦ |x|-1 ∀ k
        j
    =#
    
    
    num_sols = 0
    for i=1:tabu_list.length
        if !isnan(tabu_list.sols[i][1]) 
            num_sols += 1
        else 
            break
        end
    end

    if num_sols > 0
        M = 100
        @variable(mip_model, B[j=1:m.num_int_bin_var,k=1:num_sols], Bin)
        for k=1:num_sols, j=1:m.num_int_bin_var
            i = m.int2var_idx[j] 
            @constraint(mip_model, M*(1-B[j,k])  >= mx[i]-tabu_list.sols[k][i])
            @constraint(mip_model, -M*(1-B[j,k]) <= mx[i]-tabu_list.sols[k][i])
        end
        for k=1:num_sols
            @constraint(mip_model, sum(B[j,k] for j=1:m.num_int_bin_var) <= m.num_int_bin_var-1)
        end
        print(mip_model)
    end

    @objective(mip_model, Min, sum(mx_p[i]+mx_m[i] for i=1:m.num_int_bin_var))
    
    try 
        MathProgBase.setparameters!(m.mip_solver, TimeLimit=m.options.feasibility_pump_time_limit)
    catch
        println("Set parameters is not supported")
    end
    
    status = solve(mip_model)
    println("status: ", status)
    println("Obj: ", getobjectivevalue(mip_model))

    if num_sols > 0
        for k=1:num_sols, j=1:m.num_int_bin_var 
            println("B[",j,",",k,"]: ", getvalue(B[j,k]))
        end
        for j=1:m.num_int_bin_var
            i = m.int2var_idx[j]  
            println("mx[",i,"]: ", getvalue(mx[i]))
        end
    end

    # round mip values
    values = getvalue(mx)
    for i=1:m.num_int_bin_var
        vi = m.int2var_idx[i]
        values[vi] = round(values[vi])
    end
    return status, values
end

function generate_nlp(m, mip_sol)
    nlp_model = Model(solver=m.nl_solver)
    lb = [m.l_var; -1e6]
    ub = [m.u_var; 1e6]

    @variable(nlp_model, lb[i] <= nx[i=1:m.num_var] <= ub[i])
    setvalue(nx[1:m.num_var],mip_sol)
    # add all constraints
    for i=1:m.num_constr
        constr_expr = MathProgBase.constr_expr(m.d,i)
        expr_dereferencing!(constr_expr, nlp_model)
        JuMP.addNLconstraint(nlp_model, constr_expr)
    end

    @variable(nlp_model, nx_p[i=1:m.num_int_bin_var] >= 0)
    @variable(nlp_model, nx_m[i=1:m.num_int_bin_var] >= 0)
    for i=1:m.num_int_bin_var
        vi = m.int2var_idx[i]
        @constraint(nlp_model, nx[vi] == mip_sol[vi]+nx_p[i]-nx_m[i])
    end
    @objective(nlp_model, Min, sum(nx_p[i]+nx_m[i] for i=1:m.num_int_bin_var))

    setsolver(nlp_model, m.nl_solver)
    status = solve(nlp_model)
    println("Status: ", status)
    println("Obj: ", getobjectivevalue(nlp_model))
    return status, getvalue(nx), getobjectivevalue(nlp_model)
end


function generate_real_nlp(m,sol)
    rmodel = Model(solver=m.nl_solver)
    lb = [m.l_var; -1e6]
    ub = [m.u_var; 1e6]
    # all continuous we solve relaxation first
    @variable(rmodel, lb[i] <= rx[i=1:m.num_var] <= ub[i])
    for ni=1:m.num_int_bin_var
        i = m.int2var_idx[ni]
        JuMP.fix(rx[i], sol[i])
    end

    # define the objective function
    obj_expr = MathProgBase.obj_expr(m.d)
    expr_dereferencing!(obj_expr, rmodel)
    JuMP.setNLobjective(rmodel, m.obj_sense, obj_expr)

    # add all constraints
    for i=1:m.num_constr
        constr_expr = MathProgBase.constr_expr(m.d,i)
        expr_dereferencing!(constr_expr, rmodel)
        JuMP.addNLconstraint(rmodel, constr_expr)
    end

    status = solve(rmodel)
    return status, getvalue(rx), getobjectivevalue(rmodel)
end

function add!(t::TabuList, sol)
    t.sols[t.pointer] = sol
    t.pointer += 1
    if t.pointer > t.length
        t.pointer = 1
    end
end

function fpump(m)
    srand(1)

    start_fpump = time()
    nlp_sol = m.solution
    nlp_obj = 1 # should be not 0 for while
    c = 0
    tabu_list = TabuList()
    mip_sols = Dict{UInt64,Bool}()
    tabu_list.length = 30
    tabu_list.pointer = 1
    tabu_list.sols = []
    for i=1:tabu_list.length
        push!(tabu_list.sols, NaN*ones(length(nlp_sol))) 
    end

    fix = false
    nlp_status = :Error
    iscorrect = false
    tl = m.options.feasibility_pump_time_limit
    while !isapprox(nlp_obj,0.0, atol=atol) && time()-start_fpump < tl 
        if m.num_l_constr > 0
            mip_status, mip_sol = generate_mip(m, nlp_sol, tabu_list) 
        else
            # if no linear constraints just round the discrete variables
            mip_sol = copy(nlp_sol)
            mip_status = :Optimal
            for vi=1:m.num_int_bin_var
                vidx = m.int2var_idx[vi]
                mip_sol[vidx] = round(mip_sol[vidx])
            end
        end
        if mip_status != :Optimal
            break
        end
        add!(tabu_list, mip_sol)
        nlp_status, nlp_sol, nlp_obj = generate_nlp(m, mip_sol)
        if nlp_status != :Optimal
            break
        end
        if isapprox(nlp_obj, 0.0, atol=1e-4) && nlp_status == :Optimal
            real_status,real_sol, real_obj = generate_real_nlp(m, mip_sol)
            if real_status == :Optimal
                nlp_obj = real_obj
                nlp_sol = real_sol
                iscorrect = true
                break
            end
        end
        if haskey(mip_sols, hash(mip_sol))
            warn("Cycle detected")
        end
        for i=1:tabu_list.length
            println([tabu_list.sols[i][m.int2var_idx[j]] for j=1:m.num_int_bin_var])
        end
        mip_sols[hash(mip_sol)] = true
        c += 1
        println("c: ", c)
        if c >= 3
            error("1")
        end
    end
    
    println("It took ", time()-start_fpump, " s")
    if iscorrect
        obj = MathProgBase.eval_f(m.d, nlp_sol)
        println("obj: ", obj)
        for ni=1:m.num_int_bin_var
            i = m.int2var_idx[ni]
        end
        return nlp_sol, obj
    end
    
    return nothing, nothing
end
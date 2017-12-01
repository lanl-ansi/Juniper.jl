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

"""
    construct_affine_vector(m)

Construct a vector of affine expressions for all linear functions using the derivative
"""
function construct_affine_vector(m)
    js = MathProgBase.jac_structure(m.d)

    jg = zeros(length(js[1]))
    MathProgBase.eval_jac_g(m.d, jg, ones(m.num_var))

    # Construct the data structure for our affine constraints
    aff = Vector{Aff}(m.num_l_constr)
    for i=1:m.num_l_constr
        aff[i] = Aff()
        aff[i].var_idx = []
        aff[i].coeff = []
        constr_expr = MathProgBase.constr_expr(m.d,i)
        aff[i].rhs = constr_expr.args[3]
        aff[i].sense = constr_expr.args[1]
    end

    # if linear constraint the derivative are the coeffs
    idx = 1
    lconstr2constr = Vector{Int64}()
    constr2lconstr = Vector{Int64}()
    c = 1
    for i=1:m.num_constr
        if m.isconstrlinear[i]
            push!(lconstr2constr,i)
            push!(constr2lconstr,c)
            c+=1
        else
            push!(constr2lconstr,0)
        end
    end

    for row in js[1]
        if m.isconstrlinear[row]
            col = js[2][idx]
            aidx = constr2lconstr[row]
            push!(aff[row].var_idx, col)
            push!(aff[row].coeff, jg[idx])
        end
        idx += 1
    end
    return aff
end

"""
    generate_mip(m, nlp_sol, aff, tabu_list)

Generate a mip using the linear constraints (aff) of the original model
Minimize the distance to nlp_sol and avoid using solutions inside the tabu list
"""
function generate_mip(m, nlp_sol, aff, tabu_list)
    mip_model = Model(solver=m.mip_solver)
    lb = m.l_var
    ub = m.u_var
    @variable(mip_model, lb[i] <= mx[i=1:m.num_var] <= ub[i])
    for i=1:m.num_var
        if m.var_type[i] == :Int
            setcategory(mx[i], :Int)
        elseif m.var_type[i] == :Bin
            setcategory(mx[i], :Bin)
        end
    end
    
    for c in aff
        if c.sense == :(>=)
            @constraint(mip_model, sum(mx[c.var_idx[i]]*c.coeff[i] for i=1:length(c.var_idx)) >= c.rhs)
        elseif c.sense == :(<=)
            @constraint(mip_model, sum(mx[c.var_idx[i]]*c.coeff[i] for i=1:length(c.var_idx)) <= c.rhs)
        elseif c.sense == :(==)
            @constraint(mip_model, sum(mx[c.var_idx[i]]*c.coeff[i] for i=1:length(c.var_idx)) == c.rhs)
        end
    end

    @variable(mip_model, mabsx[i=1:m.num_int_bin_var] >= 0)
    for i=1:m.num_int_bin_var
        vi = m.int2var_idx[i]
        @constraint(mip_model, mabsx[i] >= mx[vi]-nlp_sol[vi])
        @constraint(mip_model, mabsx[i] >= -mx[vi]+nlp_sol[vi])
    end
      
    # How long is the tabu list 
    num_sols = 0
    for i=1:tabu_list.length
        if !isnan(tabu_list.sols[i][1]) 
            num_sols += 1
        else 
            break
        end
    end

    # If there solutions in the tabu list => avoid them
    if num_sols > 0
        @variable(mip_model, z1[j=1:m.num_int_bin_var,k=1:num_sols], Bin)
        @variable(mip_model, z2[j=1:m.num_int_bin_var,k=1:num_sols], Bin)
        v = tabu_list.sols
        for k=1:num_sols, j=1:m.num_int_bin_var
            i = m.int2var_idx[j]
            lbi = m.l_var[i] > typemin(Int64) ? m.l_var[i] : typemin(Int64)
            ubi = m.u_var[i] < typemax(Int64) ? m.u_var[i] : typemax(Int64)
            @constraint(mip_model, z1[j,k]+z2[j,k] <= 1)
            @constraint(mip_model, (lbi - v[k][j])*z1[j,k]+z2[j,k]+v[k][j] <= mx[i])
            @constraint(mip_model, mx[i] <= v[k][j] - z1[j,k] + (ubi-v[k][j])*z2[j,k])
        end
        for k=1:num_sols
            @constraint(mip_model, sum(z1[j,k]+z2[j,k] for j=1:m.num_int_bin_var) >= 1)
        end
    end

    @objective(mip_model, Min, sum(mabsx[i] for i=1:m.num_int_bin_var))
    
    # Break the mip solver if it takes too long or throw a warning when this option isn't available
    try 
        MathProgBase.setparameters!(m.mip_solver, TimeLimit=m.options.feasibility_pump_time_limit)
    catch
        warn("Set parameters is not supported")
    end
    
    status = solve(mip_model)
    println("MIP Status: ", status)
    println("MIP Obj: ", getobjectivevalue(mip_model))

    # round mip values
    values = getvalue(mx)
    for i=1:m.num_int_bin_var
        vi = m.int2var_idx[i]
        values[vi] = round(values[vi])
    end
    return status, values
end

"""
    generate_nlp(m, mip_sol; random_start=false)

Generates the original nlp but changes the objective to minimize the distance to the mip solution
"""
function generate_nlp(m, mip_sol; random_start=false)
    nlp_model = Model(solver=m.nl_solver)
    lb = m.l_var
    ub = m.u_var

    @variable(nlp_model, lb[i] <= nx[i=1:m.num_var] <= ub[i])
    if random_start
        for i=1:m.num_var
            lbi = m.l_var[i] > typemin(Int64) ? m.l_var[i] : typemin(Int64)
            ubi = m.u_var[i] < typemax(Int64) ? m.u_var[i] : typemax(Int64)

            if m.var_type == :Cont
                setvalue(nx[i], (ubi-lbi)*rand()+lbi)
            else
                setvalue(nx[i], rand(lbi:ubi))
            end
        end
    else
        setvalue(nx[1:m.num_var],mip_sol)
    end

    # add all constraints
    for i=1:m.num_constr
        constr_expr = MathProgBase.constr_expr(m.d,i)
        expr_dereferencing!(constr_expr, nlp_model)
        JuMP.addNLconstraint(nlp_model, constr_expr)
    end

    @variable(nlp_model, nabsx[i=1:m.num_int_bin_var] >= 0)
    for i=1:m.num_int_bin_var
        vi = m.int2var_idx[i]
        @constraint(nlp_model, nabsx[i] >= nx[vi]-mip_sol[vi])
        @constraint(nlp_model, nabsx[i] >= -nx[vi]+mip_sol[vi])
    end
    @objective(nlp_model, Min, sum(nabsx[i] for i=1:m.num_int_bin_var))

    setsolver(nlp_model, m.nl_solver)
    status = solve(nlp_model)
    println("NLP Status: ", status)
    println("NLP Obj: ", getobjectivevalue(nlp_model))
    return status, getvalue(nx), getobjectivevalue(nlp_model)
end

"""
    generate_real_nlp(m,sol)

Generate the orignal nlp and get the objective for that
"""
function generate_real_nlp(m, sol)
    rmodel = Model(solver=m.nl_solver)
    lb = m.l_var
    ub = m.u_var

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

"""
    add!(t::TabuList, m, sol)

Add a solution to the tabu list (includes only the discrete variables)
"""
function add!(t::TabuList, m, sol)
    t.sols[t.pointer] = [sol[m.int2var_idx[i]] for i=1:m.num_int_bin_var]
    t.pointer += 1
    if t.pointer > t.length
        t.pointer = 1
    end
end

"""
    fpump(m)

Run the feasibility pump 
"""
function fpump(m)
    srand(1)

    start_fpump = time()
    nlp_sol = m.solution
    nlp_obj = 1 # should be not 0 for while
    c = 1
    tabu_list = TabuList()
    mip_sols = Dict{UInt64,Bool}()
    tabu_list.length = m.options.tabu_list_length
    tabu_list.pointer = 1
    tabu_list.sols = []
    for i=1:tabu_list.length
        push!(tabu_list.sols, NaN*ones(m.num_int_bin_var)) 
    end

    # construct the linear constraints as a vector of Aff once
    if m.num_l_constr > 0
        aff = construct_affine_vector(m)
    end

    fix = false
    nlp_status = :Error
    iscorrect = false
    tl = m.options.feasibility_pump_time_limit
    while !isapprox(nlp_obj,0.0, atol=atol) && time()-start_fpump < tl 
        # generate a mip or just round if no linear constraints
        if m.num_l_constr > 0
            mip_status, mip_sol = generate_mip(m, nlp_sol, aff, tabu_list) 
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
            warn("MIP couldn't be solved to optimality")
            break
        end
        
        # If a cycle is detected which wasn't able to prevent by the tabu list (maybe too short)
        if haskey(mip_sols, hash(mip_sol))
            warn("Cycle detected")
            break
        end
        add!(tabu_list, m, mip_sol)
        mip_sols[hash(mip_sol)] = true

        nlp_status, nlp_sol, nlp_obj = generate_nlp(m, mip_sol)
        if nlp_status != :Optimal
            cnlpinf = 0 
            while cnlpinf < m.options.num_resolve_nlp_feasibility_pump && nlp_status != :Optimal
                nlp_status, nlp_sol, nlp_obj = generate_nlp(m, mip_sol; random_start=true)
                cnlpinf += 1
            end
            if nlp_status != :Optimal
                warn("NLP couldn't be solved to optimality")
                break
            end
        end
        # if the difference is near 0 => try to improve the obj by using the original obj
        if isapprox(nlp_obj, 0.0, atol=atol)
            real_status,real_sol, real_obj = generate_real_nlp(m, mip_sol)
            if real_status == :Optimal
                nlp_obj = real_obj
                nlp_sol = real_sol
                iscorrect = true
                break
            else
                nlp_obj = MathProgBase.eval_f(m.d, nlp_sol)
                iscorrect = true
                warn("Real objective wasn't solved to optimality")
            end
        end
        c += 1
    end
    
    println("It took ", time()-start_fpump, " s")
    println("It took ", c, " rounds")
    if iscorrect
        println("Obj: ", nlp_obj)
        return nlp_sol, nlp_obj
    end    
    return nothing, nothing
end
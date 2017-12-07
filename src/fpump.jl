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

function expr_dereferencing_fixing!(expr, m, var_types, sol)
    for i in 2:length(expr.args)
        if isa(expr.args[i], Union{Float64,Int64})
            k = 0
        elseif expr.args[i].head == :ref
            @assert isa(expr.args[i].args[2], Int)
            if var_types[expr.args[i].args[2]] != :Cont
                expr.args[i] = sol[expr.args[i].args[2]]
            else
                expr.args[i] = Variable(m, expr.args[i].args[2])
            end
        elseif expr.args[i].head == :call
            expr_dereferencing_fixing!(expr.args[i], m, var_types, sol)
        else
            error("expr_dereferencing :: Unexpected term in expression tree.")
        end
    end
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

    # round mip values
    values = getvalue(mx)
    for i=1:m.num_int_bin_var
        vi = m.int2var_idx[i]
        values[vi] = round(values[vi])
    end
    return status, values, getobjectivevalue(mip_model)
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
        restart_values = generate_random_restart(m)
        for i=1:m.num_var      
            setvalue(nx[i], restart_values[i])
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

    @objective(nlp_model, Min, sum((nx[m.int2var_idx[i]]-mip_sol[m.int2var_idx[i]])^2 for i=1:m.num_int_bin_var))
    setsolver(nlp_model, m.nl_solver)
    status = solve(nlp_model)
    nlp_sol = getvalue(nx)
    nx_val = getvalue(nx)
    nlp_obj = getobjectivevalue(nlp_model)
    return status, nlp_sol, nlp_obj
end

"""
    generate_real_nlp(m, sol; random_start=false)

Generate the orignal nlp and get the objective for that
"""
function generate_real_nlp(m, sol; random_start=false)
    if m.num_var == m.num_int_bin_var
        nlp_obj = MathProgBase.eval_f(m.d, sol)
        status = :Optimal
        return status, sol, nlp_obj
    end

    rmodel = Model(solver=m.nl_solver)
    lb = m.l_var
    ub = m.u_var

    @variable(rmodel, lb[i] <= rx[i=1:m.num_var] <= ub[i])
    if random_start
        restart_values = generate_random_restart(m)
        for i=1:m.num_var   
            if m.var_type[i] == :Cont   
                setvalue(rx[i], restart_values[i])
            end # discrete will be fixed anyway
        end
    end
    for i=1:m.num_int_bin_var
        vi = m.int2var_idx[i]
        JuMP.fix(rx[vi], sol[vi])
    end

    # define the objective function
    obj_expr = MathProgBase.obj_expr(m.d)
    expr_dereferencing_fixing!(obj_expr, rmodel, m.var_type, sol)
    JuMP.setNLobjective(rmodel, m.obj_sense, obj_expr)

    # add all constraints
    for i=1:m.num_constr
        constr_expr = MathProgBase.constr_expr(m.d,i)
        expr_dereferencing_fixing!(constr_expr, rmodel, m.var_type, sol)
        JuMP.addNLconstraint(rmodel, constr_expr)
    end

    status = solve(rmodel)
    real_sol = getvalue(rx)
    return status, real_sol, getobjectivevalue(rmodel)
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

function get_fp_table(mip_obj,nlp_obj,t, fields, field_chars)
    ln = ""
    i = 1
    arr = []
    for f in fields
        val = ""
        if f == :MIPobj
            val = string(round(mip_obj, 4))
        elseif f == :NLPobj
            val = string(round(nlp_obj,4))
        elseif f == :Time
            val = string(round(t,1))
        end

        if length(val) > field_chars[i]
            # too long to display shouldn't happen normally but is better than error
            # if it happens
            val = "t.l." 
        end

        padding = field_chars[i]-length(val)
        ln *= repeat(" ",trunc(Int, floor(padding/2)))
        ln *= val
        ln *= repeat(" ",trunc(Int, ceil(padding/2)))
        push!(arr,val)
        i += 1
    end
    return ln, arr
end

function print_fp_table(mip_obj,nlp_obj,t, fields, field_chars)
    ln, arr = get_fp_table(mip_obj,nlp_obj,t, fields, field_chars)
    println(ln)
end

"""
    fpump(m)

Run the feasibility pump 
"""
function fpump(m)
    srand(1)

    if are_type_correct(m.solution, m.var_type, m.int2var_idx, m.options.atol)
        return m.solution, m.objval 
    end

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

    last_table_arr = []
    fields = []
    field_chars = []
    ps = m.options.log_levels
    # Print table init
    if check_print(ps,[:Table]) 
        fields, field_chars = [:MIPobj,:NLPobj,:Time], [20,20,5]
        print_table_header(fields,field_chars)        
    end


    fix = false
    nlp_status = :Error
    iscorrect = false
    tl = m.options.feasibility_pump_time_limit
    # the tolerance can be changed => current atol
    catol = m.options.atol
    atol_counter = 0
    while !are_type_correct(nlp_sol, m.var_type, m.int2var_idx, catol) && time()-start_fpump < tl 
        # generate a mip or just round if no linear constraints
        if m.num_l_constr > 0
            mip_status, mip_sol, mip_obj = generate_mip(m, nlp_sol, aff, tabu_list) 
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
            while cnlpinf < m.options.num_resolve_nlp_feasibility_pump && nlp_status != :Optimal && time()-start_fpump < tl 
                nlp_status, nlp_sol, nlp_obj = generate_nlp(m, mip_sol; random_start=true)
                cnlpinf += 1
            end
            if nlp_status != :Optimal
                warn("NLP couldn't be solved to optimality")
                break
            end
        end

        if check_print(ps,[:Table]) 
            print_fp_table(mip_obj, nlp_obj, time()-start_fpump, fields, field_chars)
        end

        # if the current tolerance was nearly reached 5 times 
        # => If reasonable should be an option
        if atol_counter >= m.options.feasibility_pump_tolerance_counter
            catol *= 10
            warn("FPump tolerance changed to: ",catol)
            atol_counter = 0
        end

        # if the difference is near 0 => try to improve the obj by using the original obj
        # set atol for type correct to a low value as it is checked with real_nlp anyway
        if are_type_correct(nlp_sol, m.var_type, m.int2var_idx, catol*1000) || isapprox(nlp_obj, 0.0; atol=catol)
            real_status,real_sol, real_obj = generate_real_nlp(m, mip_sol)
            cnlpinf = 0
            while cnlpinf < m.options.num_resolve_nlp_feasibility_pump && real_status != :Optimal && time()-start_fpump < tl 
                real_status,real_sol, real_obj = generate_real_nlp(m, mip_sol; random_start=true)
                cnlpinf += 1
            end
            if real_status == :Optimal
                nlp_obj = real_obj
                nlp_sol = real_sol
                iscorrect = true
                break
            elseif are_type_correct(nlp_sol, m.var_type, m.int2var_idx, catol)
                nlp_obj = MathProgBase.eval_f(m.d, nlp_sol)
                iscorrect = true
                warn("Real objective wasn't solved to optimality")
                break
            end
        end
        if !isapprox(nlp_obj, 0.0; atol=catol) && isapprox(nlp_obj, 0.0; atol=10*catol)
            atol_counter += 1
        else 
            atol_counter = 0
        end
        c += 1
    end
    
    if check_print(ps,[:Table]) 
        println()
    end

    if check_print(ps,[:Info]) 
        println("FP: ", time()-start_fpump, " s")
        println("FP: ", c == 1 ? "$c round" : "$c rounds")
    end 
    m.fpump_info = Dict{Symbol,Float64}()
    m.fpump_info[:time] = time()-start_fpump
    m.fpump_info[:rounds] = c
    
    if iscorrect
        check_print(ps,[:Info]) && println("FP: Obj: ", nlp_obj)
        m.fpump_info[:obj] = nlp_obj
        m.fpump_info[:gap] = abs(m.objval-nlp_obj)/abs(nlp_obj)
        return nlp_sol, nlp_obj
    end

    m.fpump_info[:obj] = NaN
    m.fpump_info[:gap] = NaN
    check_print(ps,[:Info]) && println("FP: No integral solution found")
    return nothing, nothing
end
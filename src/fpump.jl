
type Aff
    sense     :: Symbol
    var_idx   :: Vector{Int64}
    coeff     :: Vector{Float64}
    rhs       :: Float64

    Aff() = new()
end
include("util.jl")

function generate_mip(m,nlp_sol)
    println("=========================================")
    println("generate_mip")
    println("=========================================")
    # println("sol: ", nlp_sol)
    
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
    """
        Example:
        Min x+y
        exp(x+y) <= 4
        exp(x*y) <= 4
    """

    """
        Gives the result of
        exp(x+y)-4 
        exp(x*y)-4
        evaluated at the solution
    """
    g = zeros(m.num_constr)
    MathProgBase.eval_g(m.d,g,nlp_sol)


    # println("g: ",g)

    """
        Gives the result of indices where the derivate is nonzero
        here:
        ([1,1,2,2],[1,2,1,2])
    """
    js = MathProgBase.jac_structure(m.d)

    """
        Gives the result of the derivate of g(x) at the solution
        here [10.5362, 10.5362, 4.70964, 4.70964]
        which corresponds with the indices
        ([1,1,2,2],[1,2,1,2])
        The matrix would be
        10.5 10.5
            4.7  4.7
    """
    jg = zeros(length(js[1]))
    MathProgBase.eval_jac_g(m.d, jg, nlp_sol)

    # Construct the data structure for our affine constraints
    aff = Vector{Aff}(m.num_constr)
    for i=1:m.num_constr
        aff[i] = Aff()
        aff[i].var_idx = []
        aff[i].coeff = []
        aff[i].rhs = -g[i]
        constr_expr = MathProgBase.constr_expr(m.d,i)
        aff[i].sense = constr_expr.args[1]
    end

    """
        The taylor polynom of order 1 (tangent) for exp(x+y)
        is g(sol)+g'(sol)(x-sol)
        which is g(sol)+g'(sol)x-g'(sol)*sol
    """
    idx = 1
    for row in js[1]
        col = js[2][idx]
        push!(aff[row].var_idx, col)
        push!(aff[row].coeff, jg[idx]*nlp_sol[col])
        aff[row].rhs += jg[idx]*nlp_sol[col]*nlp_sol[col]
        idx += 1
    end

    # for c in aff
        # println(c)
    # end

    
    counter = 1
    for c in aff
        if m.isconstrlinear[counter]
            # construct affine constraint 
            constr_expr = MathProgBase.constr_expr(m.d, counter)
            c = expr_linear_to_affine(constr_expr)
        end
        if c.sense == :(>=)
            @constraint(mip_model, sum(mx[c.var_idx[i]]*c.coeff[i] for i=1:length(c.var_idx)) >= c.rhs)
        elseif c.sense == :(<=)
            @constraint(mip_model, sum(mx[c.var_idx[i]]*c.coeff[i] for i=1:length(c.var_idx)) <= c.rhs)
        elseif c.sense == :(==)
            @constraint(mip_model, sum(mx[c.var_idx[i]]*c.coeff[i] for i=1:length(c.var_idx)) == c.rhs)
        else
            error(c.sense*" is not supported")
        end
        counter += 1
    end

    @variable(mip_model, mx_p[i=1:m.num_int_bin_var] >= 0)
    @variable(mip_model, mx_m[i=1:m.num_int_bin_var] >= 0)
    for i=1:m.num_int_bin_var
        vi = m.int2var_idx[i]
        @constraint(mip_model, mx[vi] == nlp_sol[vi]+mx_p[i]-mx_m[i])
    end
    @objective(mip_model, Min, sum(mx[i]-nlp_sol[m.int2var_idx[i]] for i=1:m.num_int_bin_var) 
                             + sum(nlp_sol[m.int2var_idx[i]]-mx[i] for i=1:m.num_int_bin_var)
                             + sum(mx_p[i]+mx_m[i] for i=1:m.num_int_bin_var))

    # print(mip_model)
    setsolver(mip_model, m.mip_solver)
    status = solve(mip_model)
    println("status: ", status)
    println("Obj: ", getobjectivevalue(mip_model))
    println("values: ", getvalue(mx))
    
    return getvalue(mx)
end

function generate_nlp(m, mip_sol)
    println("=========================================")
    println("generate_nlp")
    println("=========================================")
    nlp_model = Model(solver=m.nl_solver)
    lb = [m.l_var; -1e6]
    ub = [m.u_var; 1e6]

    @variable(nlp_model, lb[i] <= nx[i=1:m.num_var] <= ub[i])

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
    @objective(nlp_model, Min, sum(nx[i]-mip_sol[m.int2var_idx[i]] for i=1:m.num_int_bin_var) 
                             + sum(mip_sol[m.int2var_idx[i]]-nx[i] for i=1:m.num_int_bin_var)
                             + sum(nx_p[i]+nx_m[i] for i=1:m.num_int_bin_var))

    setsolver(nlp_model, m.nl_solver)
    status = solve(nlp_model)
    println("Status: ", status)
    println("Obj: ", getobjectivevalue(nlp_model))
    println("values: ", getvalue(nx))
    return getvalue(nx), getobjectivevalue(nlp_model)
end


function fpump(m)
    nlp_sol = m.solution
    nlp_obj = 1
    c = 0
    while !isapprox(nlp_obj,0.0, atol=atol) 
        mip_sol = generate_mip(m, nlp_sol) 
        nlp_sol, nlp_obj = generate_nlp(m, mip_sol)
        c += 1
        println("c: ", c)
    end
    
    if isapprox(nlp_obj,0.0, atol=atol)
        println("here?")
        obj = MathProgBase.eval_f(m.d, nlp_sol)
        return nlp_sol, obj
    end
    return nothing, nothing
end
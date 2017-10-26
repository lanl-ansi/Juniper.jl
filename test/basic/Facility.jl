function facility_problem()
    m = Model()
        
    file_name = "fl_10_10"
    lambda = 10
    N,M,f_s,f_c,f_pos,c_d,c_pos,dist_mat = load_fac(file_name)

    println("Init problem")
    @variable(m, sf[i=1:N], Bin)
    @variable(m, cf[c=1:M,f=1:N], Bin)
    @variable(m, md >= 0)

    @objective(m, Min, sum(sf[i]*f_s[i] for i=1:N)+lambda*md)

    @constraint(m, demand_constr[i=1:N], dot(cf[c=1:M,f=i],c_d) <= f_c[i])
    @constraint(m, served_constr[j=1:M], sum(cf[c=j,f=1:N]) == 1)
    @constraint(m, setup_constr[i=1:N], sum(cf[c=1:M,f=i]) <= M*sf[i])
    @NLconstraint(m, dist_constr[i=1:N], sum((cf[c=j,f=i]*dist_mat[j,i])^2 for j=1:M) <= md)

    return m
end
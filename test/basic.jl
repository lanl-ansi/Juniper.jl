file_name = "fl_10_10"
lambda = 10
include("load_fl.jl")

@testset "basic tests" begin


@testset "Facility" begin
    m = Model()

    println("Init problem")
    @variable(m, sf[i=1:N], Bin)
    @variable(m, cf[c=1:M,f=1:N], Bin)
    @variable(m, md >= 0)
    
    @objective(m, Min, sum(sf[i]*f_s[i] for i=1:N)+lambda*md)
    
    @constraint(m, demand_constr[i=1:N], dot(cf[c=1:M,f=i],c_d) <= f_c[i])
    @constraint(m, served_constr[j=1:M], sum(cf[c=j,f=1:N]) == 1)
    @constraint(m, setup_constr[i=1:N], sum(cf[c=1:M,f=i]) <= M*sf[i])
    @constraint(m, dist_constr[i=1:N], sum((cf[c=j,f=i]*dist_mat[j,i])^2 for j=1:M) <= md)
        
    println("Variables and constraints set")

    setsolver(m, gurobi)
    status = solve(m; relaxation=true)
    @test status == :Optimal
    println("Solved using gurobi")

    gurobi_val = getobjectivevalue(m)
    gurobi_sf = getvalue(sf)
    gurobi_cf = getvalue(cf)
    gurobi_md = getvalue(md)

    println("Set solver to MINLPBnB")
    setsolver(m, minlpbnb)
    println("Solve using MINLPBnB")
    status = solve(m)
    println("Solved using MINLPBnB")
    @test status == :Optimal

    minlpbnb_val = getobjectivevalue(m)
    minlpbnb_sf = getvalue(sf)
    minlpbnb_cf = getvalue(cf)
    minlpbnb_md = getvalue(md)
    
    println("Solution by Gurobi")
    println("obj: ", gurobi_val)
    println("sf: ", gurobi_sf)
    println("cf: ", gurobi_cf)
    println("md: ", gurobi_md)

    println("")
    println("Solution by MINLPBnb")
    println("obj: ", minlpbnb_val)
    println("sf: ", minlpbnb_sf)
    println("cf: ", minlpbnb_cf)
    println("md: ", minlpbnb_md)

    @test isapprox(minlpbnb_val, gurobi_val, atol=opt_atol, rtol=opt_rtol)
    @test isapprox(minlpbnb_sf, gurobi_sf, atol=sol_atol, rtol=sol_rtol)
    @test isapprox(minlpbnb_cf, gurobi_cf, atol=sol_atol, rtol=sol_rtol)
    @test isapprox(minlpbnb_md, gurobi_md, atol=sol_atol, rtol=sol_rtol)
end

#=@testset "PajaritoTest" begin
    m = Model(solver=minlpbnb)

    @variable(m, x >= 0, Int)
    @variable(m, y >= 0)
    @variable(m, 0 <= u <= 10, Int)
    @variable(m, w == 1)

    @objective(m, Min, -3x - y)

    @constraint(m, 3x + 10 <= 20)
    @constraint(m, y^2 <= u*w)

    status = solve(m)

    @test status == :Optimal
    @test isapprox(getobjectivevalue(m), -12.162277, atol=TOL)
    @test isapprox(getobjbound(m), -12.162277, atol=TOL)
    @test isapprox(getvalue(x), 3, atol=TOL)
    @test isapprox(getvalue(y), 3.162277, atol=TOL)
end
=#
end
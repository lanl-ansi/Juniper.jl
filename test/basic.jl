include("load_fl.jl")

@testset "basic tests" begin

#=
@testset "Facility" begin
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
        
    
    println("Variables and constraints set")

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
    
    println("")
    println("Solution by MINLPBnb")
    println("obj: ", minlpbnb_val)
    println("sf: ", minlpbnb_sf)
    println("cf: ", minlpbnb_cf)
    println("md: ", minlpbnb_md)

    @test isapprox(minlpbnb_val, 6.839268292682927, atol=opt_atol, rtol=opt_rtol)
    @test isapprox(minlpbnb_sf, [1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0], atol=sol_atol, rtol=sol_rtol)
    @test isapprox(minlpbnb_cf, [1.0 0.0 0.0 0.0 0.0 -0.0 0.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 1.0 -0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -0.0 1.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 1.0 -0.0 0.0 0.0; 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; -0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 -0.0 0.0; 0.0 0.0 0.0 0.0 1.0 -0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0]  , atol=sol_atol, rtol=sol_rtol)
    @test isapprox(minlpbnb_md, 0.08392682926829272, atol=sol_atol, rtol=sol_rtol)
end
=#


@testset "PajaritoTest" begin
    m = Model(solver=minlpbnb)

    @variable(m, x >= 0, Int)
    @variable(m, y >= 0)
    @variable(m, 0 <= u <= 10, Int)
    @variable(m, w == 1)

    @objective(m, Min, -3x - y)

    @constraint(m, 3x + 10 <= 20)
    @NLconstraint(m, y^2 <= u*w)

    status = solve(m)

    println("Obj: ", getobjectivevalue(m))
    println("x: ", getvalue(x))
    println("y: ", getvalue(y))

    @test status == :Optimal
    @test isapprox(getobjectivevalue(m), -12.162277, atol=opt_atol)
    @test isapprox(getvalue(x), 3, atol=sol_atol)
    @test isapprox(getvalue(y), 3.162277, atol=sol_atol)
end

@testset "Only one branch but several level" begin
    m = Model(solver=minlpbnb)

    @variable(m, x >= 0, Int)
    @variable(m, y >= 0, Int)
    @variable(m, 0 <= u <= 10, Int)
    @variable(m, w == 1)

    @objective(m, Min, -3x - y)

    @constraint(m, 3x + 10 <= 20)
    @NLconstraint(m, y^2 <= u*w)

    status = solve(m)

    println("Obj: ", getobjectivevalue(m))
    println("x: ", getvalue(x))
    println("y: ", getvalue(y))

    @test status == :Optimal
    @test isapprox(getobjectivevalue(m), -12, atol=opt_atol)
    @test isapprox(getvalue(x), 3, atol=sol_atol)
    @test isapprox(getvalue(y), 3, atol=sol_atol)
end


@testset "Knapsack Max" begin
    m = Model(solver=minlpbnb)

    v = [10,20,12,23,42]
    w = [12,45,12,22,21]
    @variable(m, x[1:5], Bin)

    @objective(m, Max, dot(v,x))
    
    @NLconstraint(m, sum(w[i]*x[i]^2 for i=1:5) <= 45)   

    status = solve(m)

    println("Obj: ", getobjectivevalue(m))
    println("x: ", getvalue(x))

    @test status == :Optimal
    @test isapprox(getobjectivevalue(m), 65, atol=opt_atol)
    @test isapprox(getvalue(x), [0,0,0,1,1], atol=sol_atol)
end



@testset "Facility" begin
    m = Model()

    file_name = "fl_3_4"
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
        

    setsolver(m, minlpbnb)
    status = solve(m)
    @test status == :Optimal

    minlpbnb_val = getobjectivevalue(m)
    minlpbnb_sf = getvalue(sf)
    minlpbnb_cf = getvalue(cf)
    minlpbnb_md = getvalue(md)

    println("")
    println("Solution by MINLPBnb")
    println("obj: ", minlpbnb_val)
    println("sf: ", minlpbnb_sf)
    println("cf: ", minlpbnb_cf)
    println("md: ", minlpbnb_md)

    @test isapprox(minlpbnb_val, 5.3548169043445215, atol=opt_atol, rtol=opt_rtol)
    @test isapprox(minlpbnb_sf, [1.0, 1.0, 1.0], atol=sol_atol, rtol=sol_rtol)
    @test isapprox(minlpbnb_cf, [1 0 0; 1 0 0; 0 1 0; 0 0 1] , atol=sol_atol, rtol=sol_rtol)
    @test isapprox(minlpbnb_md,  0.23548169043445205, atol=sol_atol, rtol=sol_rtol)
end


end
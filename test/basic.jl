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

#=
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

@testset "Pajarito II" begin
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
=#

@testset "Facility Small" begin
    println("==================================")
    println("FACILITY SMALL")
    println("==================================")
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



@testset "Batch.mod" begin
    println("==================================")
    println("BATCH.MOD")
    println("==================================")
    
    m = Model()
    @variable(m, 0 <= x1 <= 1.38629436111989)
    @variable(m, 0 <= x2 <= 1.38629436111989)
    @variable(m, 0 <= x3 <= 1.38629436111989)
    @variable(m, 0 <= x4 <= 1.38629436111989)
    @variable(m, 0 <= x5 <= 1.38629436111989)
    @variable(m, 0 <= x6 <= 1.38629436111989)
    @variable(m, 5.7037824746562 <= x7 <= 8.00636756765025)
    @variable(m, 5.7037824746562 <= x8 <= 8.00636756765025)
    @variable(m, 5.7037824746562 <= x9 <= 8.00636756765025)
    @variable(m, 5.7037824746562 <= x10 <= 8.00636756765025)
    @variable(m, 5.7037824746562 <= x11 <= 8.00636756765025)
    @variable(m, 5.7037824746562 <= x12 <= 8.00636756765025)
    @variable(m, 4.45966 <= x13 <= 397.747)
    @variable(m, 3.7495 <= x14 <= 882.353)
    @variable(m, 4.49144 <= x15 <= 833.333)
    @variable(m, 3.14988 <= x16 <= 638.298)
    @variable(m, 3.04452 <= x17 <= 666.667)
    @variable(m, 0.729961 <= x18 <= 2.11626)
    @variable(m, 0.530628 <= x19 <= 1.91626)
    @variable(m, 1.09024 <= x20 <= 2.47654)
    @variable(m, -0.133531 <= x21 <= 1.25276)
    @variable(m, 0.0487901 <= x22 <= 1.43508)
    @variable(m, b23, Bin)
    @variable(m, b24, Bin)
    @variable(m, b25, Bin)
    @variable(m, b26, Bin)
    @variable(m, b27, Bin)
    @variable(m, b28, Bin)
    @variable(m, b29, Bin)
    @variable(m, b30, Bin)
    @variable(m, b31, Bin)
    @variable(m, b32, Bin)
    @variable(m, b33, Bin)
    @variable(m, b34, Bin)
    @variable(m, b35, Bin)
    @variable(m, b36, Bin)
    @variable(m, b37, Bin)
    @variable(m, b38, Bin)
    @variable(m, b39, Bin)
    @variable(m, b40, Bin)
    @variable(m, b41, Bin)
    @variable(m, b42, Bin)
    @variable(m, b43, Bin)
    @variable(m, b44, Bin)
    @variable(m, b45, Bin)
    @variable(m, b46, Bin)
    
    @NLobjective(m, Min, 250*exp(x1 + 0.6*x7) + 250*exp(x2 + 0.6*x8) + 250*exp(x3 + 0.6*x9)
     + 250*exp(x4 + 0.6*x10) + 250*exp(x5 + 0.6*x11) + 250*exp(x6 + 0.6*x12))

    @NLconstraint(m, x7 - x13 >= 2.06686275947298)
    @NLconstraint(m, x8 - x13 >= 0.693147180559945)
    @NLconstraint(m, x9 - x13 >= 1.64865862558738)
    @NLconstraint(m, x10 - x13 >= 1.58923520511658)
    @NLconstraint(m, x11 - x13 >= 1.80828877117927)
    @NLconstraint(m, x12 - x13 >= 1.43508452528932)
    @NLconstraint(m, x7 - x14 >= -0.356674943938732)
    @NLconstraint(m, x8 - x14 >= -0.22314355131421)
    @NLconstraint(m, x9 - x14 >= -0.105360515657826)
    @NLconstraint(m, x10 - x14 >= 1.22377543162212)
    @NLconstraint(m, x11 - x14 >= 0.741937344729377)
    @NLconstraint(m, x12 - x14 >= 0.916290731874155)
    @NLconstraint(m, x7 - x15 >= -0.356674943938732)
    @NLconstraint(m, x8 - x15 >= 0.955511445027436)
    @NLconstraint(m, x9 - x15 >= 0.470003629245736)
    @NLconstraint(m, x10 - x15 >= 1.28093384546206)
    @NLconstraint(m, x11 - x15 >= 1.16315080980568)
    @NLconstraint(m, x12 - x15 >= 1.06471073699243)
    @NLconstraint(m, x7 - x16 >= 1.54756250871601)
    @NLconstraint(m, x8 - x16 >= 0.832909122935104)
    @NLconstraint(m, x9 - x16 >= 0.470003629245736)
    @NLconstraint(m, x10 - x16 >= 0.993251773010283)
    @NLconstraint(m, x11 - x16 >= 0.182321556793955)
    @NLconstraint(m, x12 - x16 >= 0.916290731874155)
    @NLconstraint(m, x7 - x17 >= 0.182321556793955)
    @NLconstraint(m, x8 - x17 >= 1.28093384546206)
    @NLconstraint(m, x9 - x17 >= 0.8754687373539)
    @NLconstraint(m, x10 - x17 >= 1.50407739677627)
    @NLconstraint(m, x11 - x17 >= 0.470003629245736)
    @NLconstraint(m, x12 - x17 >= 0.741937344729377)
    @NLconstraint(m, x1 + x18 >= 1.85629799036563)
    @NLconstraint(m, x2 + x18 >= 1.54756250871601)
    @NLconstraint(m, x3 + x18 >= 2.11625551480255)
    @NLconstraint(m, x4 + x18 >= 1.3609765531356)
    @NLconstraint(m, x5 + x18 >= 0.741937344729377)
    @NLconstraint(m, x6 + x18 >= 0.182321556793955)
    @NLconstraint(m, x1 + x19 >= 1.91692261218206)
    @NLconstraint(m, x2 + x19 >= 1.85629799036563)
    @NLconstraint(m, x3 + x19 >= 1.87180217690159)
    @NLconstraint(m, x4 + x19 >= 1.48160454092422)
    @NLconstraint(m, x5 + x19 >= 0.832909122935104)
    @NLconstraint(m, x6 + x19 >= 1.16315080980568)
    @NLconstraint(m, x1 + x20 >= 0)
    @NLconstraint(m, x2 + x20 >= 1.84054963339749)
    @NLconstraint(m, x3 + x20 >= 1.68639895357023)
    @NLconstraint(m, x4 + x20 >= 2.47653840011748)
    @NLconstraint(m, x5 + x20 >= 1.7404661748405)
    @NLconstraint(m, x6 + x20 >= 1.82454929205105)
    @NLconstraint(m, x1 + x21 >= 1.16315080980568)
    @NLconstraint(m, x2 + x21 >= 1.09861228866811)
    @NLconstraint(m, x3 + x21 >= 1.25276296849537)
    @NLconstraint(m, x4 + x21 >= 1.19392246847243)
    @NLconstraint(m, x5 + x21 >= 1.02961941718116)
    @NLconstraint(m, x6 + x21 >= 1.22377543162212)
    @NLconstraint(m, x1 + x22 >= 0.741937344729377)
    @NLconstraint(m, x2 + x22 >= 0.916290731874155)
    @NLconstraint(m, x3 + x22 >= 1.43508452528932)
    @NLconstraint(m, x4 + x22 >= 1.28093384546206)
    @NLconstraint(m, x5 + x22 >= 1.30833281965018)
    @NLconstraint(m, x6 + x22 >= 0.78845736036427)

    @NLconstraint(m, 250000*exp(x18 - x13) + 150000*exp(x19 - x14) + 180000*exp(x20 - x15) + 
    160000*exp(x21 - x16) + 120000*exp(x22 - x17) <= 6000)

    @NLconstraint(m,x1 - 0.693147180559945*b29 - 1.09861228866811*b35
    - 1.38629436111989*b41 == 0)

    @NLconstraint(m,x2 - 0.693147180559945*b30 - 1.09861228866811*b36
    - 1.38629436111989*b42 == 0)

    @NLconstraint(m, x3 - 0.693147180559945*b31 - 1.09861228866811*b37
    - 1.38629436111989*b43 == 0)

    @NLconstraint(m, x4 - 0.693147180559945*b32 - 1.09861228866811*b38
    - 1.38629436111989*b44 == 0)

    @NLconstraint(m, x5 - 0.693147180559945*b33 - 1.09861228866811*b39
    - 1.38629436111989*b45 == 0)

    @NLconstraint(m,  x6 - 0.693147180559945*b34 - 1.09861228866811*b40
    - 1.38629436111989*b46 == 0)

    @constraint(m,b23 + b29 + b35 + b41 == 1)
    @constraint(m,b24 + b30 + b36 + b42 == 1)
    @constraint(m,b25 + b31 + b37 + b43 == 1)
    @constraint(m,b26 + b32 + b38 + b44 == 1)
    @constraint(m,b27 + b33 + b39 + b45 == 1)
    @constraint(m,b28 + b34 + b40 + b46 == 1)

    setsolver(m, minlpbnb)
    status = solve(m)
    @test status == :Optimal

    minlpbnb_val = getobjectivevalue(m)

    println("")
    println("Solution by MINLPBnb")
    println("obj: ", minlpbnb_val)

    @test isapprox(minlpbnb_val, 285506.5082, atol=opt_atol, rtol=opt_rtol)
end


end
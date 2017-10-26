include("load_fl.jl")
include("basic/Facility.jl")
include("basic/batch.jl")

@testset "basic tests" begin

#=
@testset "Facility" begin
    m = facility_problem()
    setsolver(m, minlpbnb)
    status = solve(m)
    @test status == :Optimal

    minlpbnb_val = getobjectivevalue(m)
    minlpbnb_sf = getvalue(sf)
    minlpbnb_cf = getvalue(cf)
    minlpbnb_md = getvalue(md)

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
=#

#=
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
=#

#=
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

#=
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
=#

#=
@testset "Batch.mod" begin
    println("==================================")
    println("BATCH.MOD")
    println("==================================")
    
    m = batch_problem()

    setsolver(m, minlpbnb)
    status = solve(m)
    @test status == :Optimal

    minlpbnb_val = getobjectivevalue(m)

    println("")
    println("Solution by MINLPBnb")
    println("obj: ", minlpbnb_val)

    @test isapprox(minlpbnb_val, 285506.5082, atol=opt_atol, rtol=opt_rtol)
end
=#

#=
@testset "Batch0812_nc.mod" begin
    println("==================================")
    println("BATCH0812_NC.MOD")
    println("==================================")
    
    m = batch0812_nc_problem()

    setsolver(m, minlpbnb)
    status = solve(m)
    @test status == :Optimal

    minlpbnb_val = getobjectivevalue(m)

    println("")
    println("Solution by MINLPBnb")
    println("obj: ", minlpbnb_val)

    @test isapprox(minlpbnb_val, 2687026.6, atol=opt_atol, rtol=opt_rtol)
end
=#

#=
@testset "PowerModels case3 SOCWRPowerModel" begin
    println("==================================")
    println("case3.m")
    println("==================================")
    
    pm = build_generic_model("../PowerModels.jl/test/data/case3.m", SOCWRPowerModel, post_ots)
    m = pm.model
    @variable(m, aeiou == 1)
    @NLconstraint(m, aeiou^2 <= 1)

    setsolver(m, minlpbnb)
    status = solve(m)
    @test status == :Optimal

    minlpbnb_val = getobjectivevalue(m)

    println("")
    println("Solution by MINLPBnb")
    println("obj: ", minlpbnb_val)

    @test isapprox(minlpbnb_val, 5746.72, atol=1e0)
end
=#

#=
@testset "PowerModels case5 ACPPowerModel" begin
    println("==================================")
    println("ACPPowerModel case5.m")
    println("==================================")
    
    pm = build_generic_model("./test/data/pglib_opf_case5_pjm.m", ACPPowerModel, post_ots)
    m = pm.model
    #@variable(m, aeiou == 1)
    #@NLconstraint(m, aeiou^2 <= 1)

    setsolver(m, minlpbnb)
    status = solve(m)
    @test status == :Optimal

    minlpbnb_val = getobjectivevalue(m)

    println("")
    println("Solution by MINLPBnb")
    println("obj: ", minlpbnb_val)

    @test isapprox(minlpbnb_val, 15174.0, atol=1e0)
end
=#

#=
@testset "PowerModels case14 ACPPowerModel" begin
    println("==================================")
    println("ACPPowerModel case14.m")
    println("==================================")
    
    pm = build_generic_model("./test/data/pglib_opf_case14_ieee.m", ACPPowerModel, post_ots)
    m = pm.model

    setsolver(m, minlpbnb)
    status = solve(m)
    @test status == :Optimal

    minlpbnb_val = getobjectivevalue(m)

    println("")
    println("Solution by MINLPBnb")
    println("obj: ", minlpbnb_val)

    @test isapprox(minlpbnb_val, 6291.28, atol=1e0)
end
=#

#=
@testset "PowerModels case30 ACPPowerModel" begin
    println("==================================")
    println("ACPPowerModel case30.m")
    println("==================================")
    
    pm = build_generic_model("./test/data/pglib_opf_case30_as.m", ACPPowerModel, post_ots)
    m = pm.model

    setsolver(m, minlpbnb)
    status = solve(m)
    @test status == :Optimal

    minlpbnb_val = getobjectivevalue(m)

    println("")
    println("Solution by MINLPBnb")
    println("obj: ", minlpbnb_val)

    @test isapprox(minlpbnb_val, 8.0313e+02, atol=1e0)
end
=#

@testset "PowerModels case57 ACPPowerModel" begin
    println("==================================")
    println("ACPPowerModel case57.m")
    println("==================================")
    
    pm = build_generic_model("./test/data/pglib_opf_case57_ieee.m", ACPPowerModel, post_ots)
    m = pm.model

    setsolver(m, minlpbnb)
    status = solve(m)
    println("Status: ", status)
    @test status == :Optimal

    minlpbnb_val = getobjectivevalue(m)

    println("")
    println("Solution by MINLPBnb")
    println("obj: ", minlpbnb_val)

    @test isapprox(minlpbnb_val, 3.9323e+04, atol=1e0)
end

#=
@testset "PowerModels case6 SOCWRPowerModel" begin
    println("==================================")
    println("case6.m")
    println("==================================")
    
    pm = build_generic_model("../PowerModels.jl/test/data/case6.m", SOCWRPowerModel, post_ots)
    m = pm.model
    @variable(m, aeiou == 1)
    @NLconstraint(m, aeiou^2 <= 1)

    setsolver(m, minlpbnb)
    status = solve(m)
    @test status == :Optimal

    minlpbnb_val = getobjectivevalue(m)

    println("")
    println("Solution by MINLPBnb")
    println("obj: ", minlpbnb_val)

    @test isapprox(minlpbnb_val, 11559.8, atol=1e0)
end
=#

#=
@testset "PowerModels case5 SOCWRPowerModel" begin
    println("==================================")
    println("case5.m")
    println("==================================")
    
    pm = build_generic_model("../PowerModels.jl/test/data/case5.m", SOCWRPowerModel, post_ots)
    m = pm.model
    @variable(m, aeiou == 1)
    @NLconstraint(m, aeiou^2 <= 1)

    setsolver(m, minlpbnb)
    status = solve(m)
    @test status == :Optimal

    minlpbnb_val = getobjectivevalue(m)

    println("")
    println("Solution by MINLPBnb")
    println("obj: ", minlpbnb_val)

    @test isapprox(minlpbnb_val, 14999.7, atol=1e0)
end
=#

#=
@testset "PowerModels case14 SOCWRPowerModel" begin
    println("==================================")
    println("case14.m")
    println("==================================")
    
    pm = build_generic_model("../PowerModels.jl/test/data/case14.m", ACPPowerModel, post_ots)
    m = pm.model
    @variable(m, aeiou == 1)
    @NLconstraint(m, aeiou^2 <= 1)

    setsolver(m, minlpbnb)
    status = solve(m)
    @test status == :Optimal

    minlpbnb_val = getobjectivevalue(m)

    println("")
    println("Solution by MINLPBnb")
    println("obj: ", minlpbnb_val)

    @test isapprox(minlpbnb_val, 8075.1, atol=1e0)
end
=#

end
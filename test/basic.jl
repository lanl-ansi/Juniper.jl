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
    @NLconstraint(m, dist_constr[i=1:N], sum((cf[c=j,f=i]*dist_mat[j,i])^2 for j=1:M) <= md)
        
    
    println("Variables and constraints set")

    #=
    setsolver(m, gurobi)
    status = solve(m; relaxation=true)

    gurobi_val = getobjectivevalue(m)
    gurobi_sf = getvalue(sf)
    gurobi_cf = getvalue(cf)
    gurobi_md = getvalue(md)

    println("Solution by Gurobi")
    println("obj: ", gurobi_val)
    println("sf: ", gurobi_sf)
    println("cf: ", gurobi_cf)
    println("md: ", gurobi_md)

    @test status == :Optimal
    println("Solved using gurobi")
    =#
 

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

    @test isapprox(minlpbnb_val, 1.1032368700581892, atol=opt_atol, rtol=opt_rtol)
    @test isapprox(minlpbnb_sf, [0.159463, 0.0923787, 0.115921, 0.103829, 0.0821798, 0.0882018, 0.0846488, 0.142793, 0.0773032, 0.0532822], atol=sol_atol, rtol=sol_rtol)
    @test isapprox(minlpbnb_cf, [0.559328 0.172821 0.0701434 0.0367462 0.0241976 0.031867 0.0311984 0.0255741 0.0252636 0.0228609; 0.653452 0.179854 0.0594497 0.0283198 0.0176734 0.0119706 0.0128203 0.0116685 0.012603 0.0121886; 0.207836 0.190494 0.105082 0.0607358 0.0415738 0.109111 0.0990059 0.0711929 0.0627655 0.0522035; 0.0149891 0.0300479 0.0352636 0.0319615 0.0265611 0.0554079 0.157491 0.329317 0.215524 0.103437; 0.0125256 0.0185618 0.0175857 0.0151022 0.0135328 0.44569 0.253081 0.105814 0.069241 0.0488658; 0.023643 0.106172 0.543573 0.212845 0.0618376 0.00683552 0.00892605 0.0099389 0.0126688 0.0135605; 0.0829059 0.118753 0.086607 0.0555725 0.0392962 0.174572 0.178026 0.113428 0.0863921 0.0644469; 0.0187267 0.0462243 0.0870137 0.173171 0.314166 0.0216304 0.0353812 0.0544523 0.101425 0.147809; 0.00238136 0.00431537 0.00489697 0.00480892 0.00457633 0.0161185 0.0586521 0.692395 0.16749 0.0443656; 0.0188393 0.0565448 0.149593 0.419023 0.278383 0.00881397 0.0119066 0.0141529 0.0196595 0.0230841], atol=sol_atol, rtol=sol_rtol)
    @test isapprox(minlpbnb_md, 0.01032368700571742, atol=sol_atol, rtol=sol_rtol)
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
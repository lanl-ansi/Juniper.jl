
@testset "Function testing" begin

@testset ":Options" begin
    m = Model(solver=DefaultTestSolver(;traverse_strategy=:DBFS,obj_epsilon=0.5))

    v = [10,20,12,23,42]
    w = [12,45,12,22,21]
    @variable(m, x[1:5], Bin)

    @objective(m, Max, dot(v,x))

    @NLconstraint(m, sum(w[i]*x[i]^2 for i=1:5) <= 45)   

    JuMP.build(m)
    options = m.internalModel.options

    nd_options = MINLPBnB.get_non_default_options(options)    
    @test nd_options[:obj_epsilon] == 0.5
    @test nd_options[:traverse_strategy] == :DBFS
    @test length(nd_options[:log_levels]) == 0
end

@testset "Table config" begin
    m = Model(solver=DefaultTestSolver(;branch_strategy=:StrongPseudoCost,processors=2,
                                        strong_restart=true))

    v = [10,20,12,23,42]
    w = [12,45,12,22,21]
    @variable(m, x[1:5], Bin)

    @objective(m, Max, dot(v,x))

    @NLconstraint(m, sum(w[i]*x[i]^2 for i=1:5) <= 45)   
    
    JuMP.build(m)
    m = m.internalModel
    options = m.options

    fields, field_chars = MINLPBnB.get_table_config(options)
    start_time = time()
    tree = MINLPBnB.init(start_time,m)
    node = MINLPBnB.new_default_node(1,1,zeros(Int64,5),ones(Int64,5),zeros(Int64,5))
    step_obj = MINLPBnB.new_default_step_obj(m,node)
    counter = 1
    tab_ln, tab_arr = MINLPBnB.get_table_line(2,tree,node,step_obj,start_time,fields,field_chars,counter;last_arr=[])
    @test length(fields) == length(field_chars)
    for i in 1:length(fields)
        fc = field_chars[i]
        f = string(fields[i])
        @test fc >= length(f)
    end
    @test length(fields) == length(tab_arr)
    @test length(tab_ln) == sum(field_chars)
end

end
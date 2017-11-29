
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

    nd_options = Juniper.get_non_default_options(options)    
    @test nd_options[:obj_epsilon] == 0.5
    @test nd_options[:traverse_strategy] == :DBFS
    @test length(nd_options[:log_levels]) == 0
end

function option_not_available_t()
    m = Model(solver=DefaultTestSolver(;traverse_strategy=:DBS,obj_epsilon=0.5))
end
function option_not_available_b()
    m = Model(solver=DefaultTestSolver(;branch_strategy=:Pseudo,obj_epsilon=0.5))
end
function option_not_available()
    m = Model(solver=DefaultTestSolver(;branch=:Pseudo,obj_epsilon=0.5))
end


@testset ":Option not available" begin
    @test_throws ErrorException option_not_available_t()
    @test_throws ErrorException option_not_available_b()
    @test !isa(try option_not_available() catch ex ex end, Exception) 
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

    fields, field_chars = Juniper.get_table_config(options)
    ln = Juniper.get_table_header_line(fields, field_chars)
    @test length(ln) == sum(field_chars)
        
    start_time = time()
    tree = Juniper.init(start_time,m)
    tree.best_bound = 100
    node = Juniper.new_default_node(1,1,zeros(Int64,5),ones(Int64,5),zeros(Int64,5))
    step_obj = Juniper.new_default_step_obj(m,node)
    step_obj.counter = 1
    tab_ln, tab_arr = Juniper.get_table_line(2,tree,node,step_obj,start_time,fields,field_chars;last_arr=[])
    @test length(fields) == length(field_chars)
    for i in 1:length(fields)
        fc = field_chars[i]
        f = string(fields[i])
        @test fc >= length(f)
    end
    # Test with incumbent
    tree.incumbent = Juniper.Incumbent(42,[0,0,0,0,1],:UserLimit,65)
    tab_ln, tab_arr = Juniper.get_table_line(2,tree,node,step_obj,start_time,fields,field_chars;last_arr=[])
    @test length(fields) == length(field_chars)
    for i in 1:length(fields)
        fc = field_chars[i]
        f = string(fields[i])
        @test fc >= length(f)
    end
    @test length(fields) == length(tab_arr)
    @test length(tab_ln) == sum(field_chars)

    # test for diff
    tab_arr_new = copy(tab_arr)
    tab_arr_new[1] = 100
    @test Juniper.is_table_diff(fields, tab_arr, tab_arr_new) == true
    @test Juniper.is_table_diff(fields, tab_arr, tab_arr) == false
    
    # test if :Time not exists
    i = 1
    for f in fields
        if f == :Time
            deleteat!(fields,i)
            deleteat!(field_chars,i)
        end
        i+= 1
    end
    tab_ln, tab_arr = Juniper.get_table_line(2,tree,node,step_obj,start_time,fields,field_chars;last_arr=[])
    tab_arr_new = copy(tab_arr)
    tab_arr_new[1] = 100
    @test Juniper.is_table_diff(field_chars, tab_arr, tab_arr_new) == true
    @test Juniper.is_table_diff(field_chars, tab_arr, tab_arr) == false

    # produce table line with high counter for nrestarts
    step_obj.counter = 100
    tab_ln, tab_arr = Juniper.get_table_line(2,tree,node,step_obj,start_time,fields,field_chars;last_arr=[])
    idx_restarts = findfirst(fields .== :Restarts)
    @test tab_arr[idx_restarts] == "-"
    
    # normal gain gap
    step_obj.gain_gap = 0.05
    tab_ln, tab_arr = Juniper.get_table_line(2,tree,node,step_obj,start_time,fields,field_chars;last_arr=[])
    idx_gain_gap = findfirst(fields .== :GainGap)
    @test tab_arr[idx_gain_gap] == "5.0%" || tab_arr[idx_gain_gap] == "5.00%"

    # Inf gain gap
    step_obj.gain_gap = Inf
    tab_ln, tab_arr = Juniper.get_table_line(2,tree,node,step_obj,start_time,fields,field_chars;last_arr=[])
    idx_gain_gap = findfirst(fields .== :GainGap)
    @test tab_arr[idx_gain_gap] == "∞"

    # Huge gain gap
    step_obj.gain_gap = 123456789000
    tab_ln, tab_arr = Juniper.get_table_line(2,tree,node,step_obj,start_time,fields,field_chars;last_arr=[])
    idx_gain_gap = findfirst(fields .== :GainGap)
    @test tab_arr[idx_gain_gap] == ">>"

end

end
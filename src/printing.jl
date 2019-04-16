"""
    check_print(vec::Vector{Symbol}, ps::Vector{Symbol})

Check whether something should get printed. Only if the log_levels option
include the necessary level 
"""
function check_print(log_levels::Vector{Symbol}, necessary::Vector{Symbol})
    for v in log_levels
        if v in necessary
            return true
        end
    end
    return false
end

function print_info(m::JuniperProblem)
    println("#Variables: ", m.num_var)
    println("#IntBinVar: ", m.num_disc_var)
    println("#Constraints: ", m.num_constr)
    println("#Linear Constraints: ", m.num_l_constr)
    println("#Quadratic Constraints: ", m.num_q_constr)
    println("#NonLinear Constraints: ", m.num_nl_constr)
    println("Obj Sense: ", m.obj_sense)
    println()
end

function print_dict(d)
    longest_key_name = maximum([length(string(key)) for key in keys(d)])+2
    for key in keys(d)
        skey = string(key)
        pkey = skey*repeat(" ", longest_key_name-length(skey))
        if d[key] == nothing
            println(pkey, ": NA")
        else
            println(pkey, ": ",d[key])
        end
    end
end

function get_non_default_options(options)
    defaults = Juniper.get_default_options()
    non_defaults = Dict{Symbol,Any}()
    for fname in fieldnames(SolverOptions)
        # TODO: Better printing of nl_solver/mip_solver name when SolverName exists

        # doesn't work for arrays but the only array atm is log_levels 
        # and the default doesn't include :Options therefore !== should work...
        if getfield(options,fname) !== getfield(defaults,fname)
            non_defaults[fname] = getfield(options,fname)
        end
    end
    return non_defaults
end

function print_options(m::JuniperProblem; all=true)
    if all
        println(m.options)
    else
        print_dict(get_non_default_options(m.options))
    end
    println()
end

function print_fp_table(mip_obj,nlp_obj,t, fields, field_chars, catol)
    ln, arr = get_fp_table(mip_obj,nlp_obj,t, fields, field_chars, catol)
    println(ln)
end

function print_final_timing(time_bnb_solve::Float64, time_obj::TimeObj)
    println("BnB time: ", round(time_bnb_solve; digits=2))
    println("% solve child time: ", round((time_obj.solve_leaves_get_idx+time_obj.solve_leaves_branch)/time_bnb_solve*100; digits=1))
    println("Solve node time get idx: ", round(time_obj.solve_leaves_get_idx; digits=2))
    println("Solve node time branch: ", round(time_obj.solve_leaves_branch; digits=2))
    println("Branch time: ", round(time_obj.branch; digits=2))
    println("Get idx time: ", round(time_obj.get_idx; digits=2))
    println("Upd gains time: ", round(time_obj.upd_gains; digits=2))
end
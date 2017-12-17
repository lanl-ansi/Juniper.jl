function get_table_config(opts)
    fields = [:ONodes,:CLevel,:Incumbent,:BestBound,:Gap,:Time]
    field_chars = [9,8,28,28,7,8]
    opbs = opts.branch_strategy

    if (opbs == :StrongPseudoCost || opbs == :Reliability) && opts.strong_restart
        push!(fields,:Restarts)
        push!(field_chars,10)
    end

    if opts.processors > 1
        unshift!(fields,:p)
        unshift!(field_chars, 3)
    end

    if opbs == :StrongPseudoCost || opbs == :PseudoCost || opbs == :Reliability
        push!(fields, :GainGap)
        push!(field_chars, 10)
    end
    return fields, field_chars
end

function print_table_header(fields, field_chars)
    println("")
    println(get_table_header_line(fields, field_chars))
    println(repeat("=", sum(field_chars)))
end

function get_table_header_line(fields, field_chars) 
    ln = ""
    i = 1
    for f in fields
        f = string(f)
        padding = field_chars[i]-length(f)
        ln *= repeat(" ",trunc(Int, floor(padding/2)))
        ln *= f
        ln *= repeat(" ",trunc(Int, ceil(padding/2)))
        i += 1
    end
    return ln
end

function is_table_diff(fields, last_arr, new_arr)
    if length(last_arr) != length(new_arr)
        return true
    end    

    time_idx = findfirst(fields .== :Time)
    if time_idx != 0
        last_arr = vcat(last_arr[1:time_idx-1],last_arr[time_idx+1:end])
        new_arr = vcat(new_arr[1:time_idx-1],new_arr[time_idx+1:end])
    end
    for i=1:length(last_arr)
        last_arr[i] != new_arr[i] && return true
    end
    return false 
end

function print_table(p, tree, node, step_obj, start_time, fields, field_chars; last_arr=[])
    table_line, table_arr = get_table_line(p, tree, node, step_obj, start_time, fields, field_chars;
                                           last_arr=[])    
    is_table_diff(fields, last_arr,table_arr) && println(table_line)
    return table_arr
end

function get_table_line(p, tree, node, step_obj, start_time, fields, field_chars; last_arr=[])
    gain_gap = step_obj.gain_gap
    opbs = tree.options.branch_strategy
    counter = step_obj.counter
    if (opbs != :StrongPseudoCost || counter > tree.options.strong_branching_nsteps) && step_obj.upd_gains != :GainsToTree
        step_obj.nrestarts = -1 # will be displayed as -
    end
    if counter <= tree.options.strong_branching_nsteps || step_obj.upd_gains == :GainsToTree
        gain_gap = -1.0 # will be displayed as -
    end

    arr = []
    nrestarts = step_obj.nrestarts
    i = 1
    ln = ""
    for f in fields
        val = ""
        if f == :p
            val = string(p)
        elseif f == :ONodes
            val = string(length(tree.branch_nodes))
        elseif f == :Incumbent
            if isdefined(tree,:incumbent) 
                incu = tree.incumbent
                val = string(round(incu.objval,2))
            else 
                val = "-"
            end
        elseif f == :BestBound
            val = string(round(tree.best_bound,2))
        elseif f == :Gap
            if isdefined(tree,:incumbent)
                b = tree.best_bound
                o = tree.incumbent.objval
                val = round(abs(b-o)/abs(o)*100,2)
                if length(string(val)) > field_chars[i]
                    if val > 0 && val < tree.options.mip_gap*100
                        val = "< "*string(tree.options.mip_gap*100)*"%"
                    elseif val > 1000
                        val = ">>"
                    else
                        val = string(val)*"%"
                    end
                else
                    val = string(val)*"%"
                end
            else
                val = "-"
            end
        elseif f == :Time
            val = string(round(time()-start_time,1))
        elseif f == :Restarts
            if nrestarts == -1
                val = "-"
            else
                val = string(nrestarts)
            end
        elseif f == :CLevel
            val = string(node.level+1)
        elseif f == :GainGap
            if gain_gap == Inf
                val = "∞"
            elseif gain_gap == -1.0
                val = "-"
            else
                val = string(round(gain_gap*100,1))*"%"
            end
            if length(val) > field_chars[i]
                val = ">>"
            end
        end
        if length(val) > field_chars[i]
            # too long to display shouldn't happen normally but is better than error
            # if it happens
            val = "t.l." 
        end

        padding = field_chars[i]-length(val)
        ln *= repeat(" ",trunc(Int, floor(padding/2)))
        ln *= val
        ln *= repeat(" ",trunc(Int, ceil(padding/2)))
        push!(arr,val)
        i += 1
    end
    return ln, arr
end
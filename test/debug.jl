#=
    Provides different functions to check the debugTree
=#

"""
    traverse(entry,callback,params,results)

Traverse over the debug dict and call callback for every entry with the hash, step_obj
as well as parameters.
Returns the results array in the end.
"""
function traverse(entry,callback,params,results)
    push!(results,callback(entry[:hash],entry[:step_obj],params))
    if haskey(entry,:children)
        results = traverse(entry[:children][1], callback, params, results)
        results = traverse(entry[:children][2], callback, params, results)
    end
    return results
end

"""
    c_isstate(hash,entry,params)

A callback function for traverse which checks whether a given state (params[:state]) is present
"""
function c_isstate(hash,entry,params)
    return entry[:node][:state] == params[:state]
end

"""
    c_hashes(hash,entry,params)

A callback function for traverse which gets the hash
"""
function c_hashes(hash,entry,params)
    return hash
end

"""
    c_counter(hash,entry,params)

A callback function for traverse which gets the counter
"""
function c_counter(hash,entry,params)
    return haskey(entry,:counter) ? entry[:counter] : 0
end


"""
    getnstate(d,state)

Get number of how often a given state arised
"""
function getnstate(d,state)
    dictTree = d[:tree]
    params = Dict{Symbol,Symbol}()
    params[:state] = state
    results = Vector{Bool}()
    results = traverse(dictTree, c_isstate, params, results)
    return sum(results)
end

"""
    different_hashes(d)

Return if all hashes are different
"""
function different_hashes(d)
    dictTree = d[:tree]
    results = Vector{String}()
    results = traverse(dictTree, c_hashes, [], results)
    return length(Set(results)) == length(results)
end

"""
    counter_test(d,nbranches)

Check whether every counter is used exactly onces and that it is the same as the 
number of branches.
"""
function counter_test(d, nbranches)
    dictTree = d[:tree]
    results = Vector{Int64}()
    results = traverse(dictTree, c_counter, [], results)
    sort!(results)
    results = results[something(findlast(isequal(0), results), 0)+1:end]
    @assert length(Set(results)) == length(results)
    @assert results[1] == 1
    @assert results[end] == length(results)
end

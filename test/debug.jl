"""
    Provides different functions to check the debugTree
"""

function traverse(node,callback,params,results)
    push!(results,callback(node[:node],params))
    if haskey(node,:children)
        results = traverse(node[:children][1], callback, params, results)
        results = traverse(node[:children][2], callback, params, results)
    end
    return results
end

function isState(node,params)
    return node[:state] == params[:state]
end

function getnState(d,state)
    dictTree = d[:tree]
    params = Dict{Symbol,Symbol}()
    params[:state] = state
    results = Vector{Bool}()
    results = traverse(dictTree, isState, params, results)
    return sum(results)
end
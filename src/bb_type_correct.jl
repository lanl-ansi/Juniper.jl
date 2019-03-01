"""
    is_type_correct(x, var_type, atol)

Check whether a variable x has the correct type
"""
function is_type_correct(x, var_type, atol)
    if var_type != :Cont
        if !isapprox(round(x)-x, 0; atol=atol)
            return false
        end
    end
    return true
end

"""
    are_type_correct(sol, types, disc2var_idx, atol)

Check whether all variables have the correct type
"""
function are_type_correct(sol, types, disc2var_idx, atol)
    for i in disc2var_idx
        if !isapprox(round(sol[i])-sol[i], 0; atol=atol)
            return false
        end
    end
    return true
end

"""
    all_reasonable_type_correct(sol, disc2var_idx, reasonable_idx, atol)

Check whether all reasonable variable have are discrete alreadyy
"""
function all_reasonable_type_correct(sol, disc2var_idx, reasonable_idx, atol)
    for i in reasonable_idx
        idx = disc2var_idx[i]
        if !isapprox(round(sol[idx])-sol[idx], 0; atol=atol)
            return false
        end
    end
    return true
end
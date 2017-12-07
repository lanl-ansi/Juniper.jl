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
are_type_correct(sol, types, int2var_idx, atol)

Check whether all variables have the correct type
"""
function are_type_correct(sol, types, int2var_idx, atol)
    for i in int2var_idx
        if !isapprox(round(sol[i])-sol[i], 0; atol=atol)
            return false
        end
    end
    return true
end
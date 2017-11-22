"""
    is_type_correct(x, var_type)

Check whether a variable x has the correct type
"""
function is_type_correct(x, var_type)
    if var_type != :Cont
        if !isapprox(abs(round(x)-x),0, atol=atol, rtol=rtol)
        return false
        end
    end
    return true
end

"""
are_type_correct(sol, types, int2var_idx)

Check whether all variables have the correct type
"""
function are_type_correct(sol, types, int2var_idx)
    for i in int2var_idx
        if !isapprox(round(sol[i])-sol[i],0, atol=atol, rtol=rtol)
            return false
        end
    end
    return true
end
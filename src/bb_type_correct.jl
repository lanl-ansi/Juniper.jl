"""
    is_type_correct(x,var_type)

Check whether a variable x has the correct type
"""
function is_type_correct(x,var_type)
    if var_type != :Cont
        if !isapprox(abs(round(x)-x),0, atol=atol, rtol=rtol)
        return false
        end
    end
    return true
end

"""
are_type_correct(x,types)

Check whether all variables have the correct type
"""
function are_type_correct(x,types)
    for i=1:length(x)
        if types[i] != :Cont
            if !isapprox(abs(round(x[i])-x[i]),0, atol=atol, rtol=rtol)
                return false
            end
        end
    end
    return true
end
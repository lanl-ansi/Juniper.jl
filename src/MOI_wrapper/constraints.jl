"""
MOI constraints
"""

"""
Single variable bound constraints 
"""
MOI.supports_constraint(::Optimizer, ::Type{SVF}, ::Type{MOI.LessThan{Float64}}) = true
MOI.supports_constraint(::Optimizer, ::Type{SVF}, ::Type{MOI.GreaterThan{Float64}}) = true
MOI.supports_constraint(::Optimizer, ::Type{SVF}, ::Type{MOI.EqualTo{Float64}}) = true

"""
Binary/Integer variable support
"""
MOI.supports_constraint(::Optimizer, ::Type{SVF}, ::Type{<:VAR_TYPES}) = true

"""
Linear constraints
"""
MOI.supports_constraint(::Optimizer, ::Type{SAF}, ::Type{<:BOUNDS}) = true

"""
Quadratic constraints (scalar i.e., vectorized constraints are not supported)
"""
MOI.supports_constraint(::Optimizer, ::Type{SQF}, ::Type{MOI.LessThan{Float64}}) = true
MOI.supports_constraint(::Optimizer, ::Type{SQF}, ::Type{MOI.GreaterThan{Float64}}) = true
MOI.supports_constraint(::Optimizer, ::Type{SQF}, ::Type{MOI.EqualTo{Float64}}) = true


"""
``define_add_constraint`` macro - from Ipopt.jl
"""
macro define_add_constraint(function_type, set_type, array_name)
    quote
        function MOI.add_constraint(model::Optimizer, func::$function_type, set::$set_type)
            check_inbounds(model, func)
            push!(model.$(array_name), (func, set))
            return MOI.ConstraintIndex{$function_type, $set_type}(length(model.$(array_name)))
        end
    end
end

"""
``MOI.add_constraint()`` overloads for all the supported constraint types 
"""
@define_add_constraint(SAF, MOI.LessThan{Float64}, linear_le_constraints)
@define_add_constraint(SAF, MOI.GreaterThan{Float64}, linear_ge_constraints)
@define_add_constraint(SAF, MOI.EqualTo{Float64}, linear_eq_constraints)
@define_add_constraint(SQF, MOI.LessThan{Float64}, quadratic_le_constraints)
@define_add_constraint(SQF, MOI.GreaterThan{Float64}, quadratic_ge_constraints)
@define_add_constraint(SQF, MOI.EqualTo{Float64}, quadratic_eq_constraints)

"""
Binary variable support 
"""
function MOI.add_constraint(model::Optimizer, v::SVF, ::MOI.ZeroOne)
    vi = v.variable
    check_inbounds(model, vi)
    # the bounds are set using info_array_of_variables 
    # according to mlubin the bounds should not be set here as the bounds should stay when 
    # the binary constraint is deleted
    model.variable_info[vi.value].is_binary = true
	
    return MOI.ConstraintIndex{SVF, MOI.ZeroOne}(vi.value)
end 

"""
Integer variable support 
"""
function MOI.add_constraint(model::Optimizer, v::SVF, ::MOI.Integer)
    vi = v.variable
	check_inbounds(model, vi)
    model.variable_info[vi.value].is_integer = true
	
    return MOI.ConstraintIndex{SVF, MOI.Integer}(vi.value)
end 

"""
ConstraintIndex support 
"""
MOI.supports(::Optimizer, ::MOI.ConstraintName, ::Type{CI}) = true
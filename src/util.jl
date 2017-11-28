"""
This function takes a constraint/objective expression and converts it into a affine expression data structure
Use the function to traverse linear expressions traverse_expr_linear_to_affine()
"""
function expr_linear_to_affine(expr)
    # The input should follow :(<=, LHS, RHS)
    aff = Aff()
    if expr.args[1] in [:(==), :(>=), :(<=)] # For a constraint expression
        @assert isa(expr.args[3], Float64) || isa(expr.args[3], Int)
        @assert isa(expr.args[2], Expr)
        # non are buffer spaces, not used anywhere
        lhscoeff, lhsvars, rhs, non, non = traverse_expr_linear_to_affine(expr.args[2])
        rhs = -rhs + expr.args[3]
        aff.sense = expr.args[1]
    elseif expr.head == :ref  # For single variable objective expression
        lhscoeff = [1.0]
        lhsvars = [expr.args[2]]
        rhs = 0
        aff.sense = nothing
    else # For an objective expression
        lhscoeff, lhsvars, rhs, non, non = traverse_expr_linear_to_affine(expr)
        aff.sense = nothing
    end

    
    aff.coeff = lhscoeff
    aff.var_idx	= lhsvars
    aff.rhs	= rhs

    return aff
end

"""
traverse_expr_linear_to_affine(expr, lhscoeffs=[], lhsvars=[], rhs=0.0, bufferVal=0.0, bufferVar=nothing, sign=1.0, level=0)
This function traverse a left hand side tree to collect affine terms.
Updated status : possible to handle (x-(x+y(t-z))) cases where signs are handled properly
"""
function traverse_expr_linear_to_affine(expr, lhscoeffs=[], lhsvars=[], rhs=0.0, bufferVal=0.0, bufferVar=nothing, sign=1.0, coef=1.0, level=0)
    # @show expr, coef, bufferVal

    reversor = Dict(true => -1.0, false => 1.0)
    function sign_convertor(subexpr, pos)
        if length(subexpr.args) == 2 && subexpr.args[1] == :-
            return -1.0
        elseif length(subexpr.args) > 2 && subexpr.args[1] == :- && pos > 2
            return -1.0
        end
        return 1.0
    end

    if isa(expr, Float64) || isa(expr, Int) # Capture any coefficients or right hand side
        (bufferVal > 0.0) ? bufferVal *= expr : bufferVal = expr * coef
        return lhscoeffs, lhsvars, rhs, bufferVal, bufferVar
    elseif expr in [:+, :-]
        if bufferVal != 0.0 && bufferVar != nothing
            push!(lhscoeffs, bufferVal)
            push!(lhsvars, bufferVar)
            bufferVal = 0.0
            bufferVar = nothing
        end
        return lhscoeffs, lhsvars, rhs, bufferVal, bufferVar
    elseif expr in [:*]
        return lhscoeffs, lhsvars, rhs, bufferVal, bufferVar
    elseif expr in [:(<=), :(==), :(>=)]
        return lhscoeffs, lhsvars, rhs, bufferVal, bufferVar
    elseif expr in [:/, :^]
        error("Unsupported operators $expr, it is suppose to be affine function")
    elseif expr.head == :ref
        bufferVar = expr.args[2]
        return lhscoeffs, lhsvars, rhs, bufferVal, bufferVar
    end

    # HOTPATCH : Special Structure Recognization
    start_pos = 1
    if (expr.args[1] == :*) && (length(expr.args) == 3)
        if (isa(expr.args[2], Float64) || isa(expr.args[2], Int)) && (expr.args[3].head == :call)
            (coef != 0.0) ? coef = expr.args[2] * coef : coef = expr.args[2]	# Patch
            # coef = expr.args[2]
            start_pos = 3
            warn("Speicial expression structure detected [*, coef, :call, ...]. Currently handling using a beta fix...", once=true)
        end
    end

    for i in start_pos:length(expr.args)
        lhscoeff, lhsvars, rhs, bufferVal, bufferVar = traverse_expr_linear_to_affine(expr.args[i], lhscoeffs, lhsvars, rhs, bufferVal, bufferVar, sign*sign_convertor(expr, i), coef, level+1)
        if expr.args[1] in [:+, :-]  # Term segmentation [:-, :+], see this and wrap-up the current (linear) term
            if bufferVal != 0.0 && bufferVar != nothing  # (sign) * (coef) * (var) => linear term
                push!(lhscoeffs, sign*sign_convertor(expr, i)*bufferVal)
                push!(lhsvars, bufferVar)
                bufferVal = 0.0
                bufferVar = nothing
            end
            if bufferVal != 0.0 && bufferVar == nothing  # (sign) * (coef) => right-hand-side term
                rhs += sign*sign_convertor(expr, i)*bufferVal
                bufferVal = 0.0
            end
            if bufferVal == 0.0 && bufferVar != nothing && expr.args[1] == :+
                push!(lhscoeffs, sign*1.0*coef)
                push!(lhsvars, bufferVar)
                bufferVar = nothing
            end
            if bufferVal == 0.0 && bufferVar != nothing && expr.args[1] == :-
                push!(lhscoeffs, sign*sign_convertor(expr, i)*coef)
                push!(lhsvars, bufferVar)
                bufferVar = nothing
            end
        elseif expr.args[1] in [:(<=), :(==), :(>=)]
            rhs = expr.args[end]
        end
    end

    if level == 0
        if bufferVal != 0.0 && bufferVar != nothing
            push!(lhscoeffs, bufferVal)
            push!(lhsvars, bufferVar)
            bufferVal = 0.0
            bufferVar = nothing
        end
    end

    return lhscoeffs, lhsvars, rhs, bufferVal, bufferVar
end
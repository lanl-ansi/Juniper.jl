abstract type AbstractModelFilter <: MOI.ModelLike end

function MOI.get(lf::AbstractModelFilter, attr::MOI.AbstractModelAttribute)
    return MOI.get(lf.inner, attr)
end
function MOI.get(
    lf::AbstractModelFilter,
    attr::MOI.AbstractVariableAttribute,
    vi::MOI.VariableIndex,
)
    return MOI.get(lf.inner, attr, vi)
end
function MOI.get(
    lf::AbstractModelFilter,
    attr::MOI.AbstractConstraintAttribute,
    ci::MOI.ConstraintIndex,
)
    return MOI.get(lf.inner, attr, ci)
end

struct NoObjectiveFilter{M<:MOI.ModelLike} <: AbstractModelFilter
    inner::M
end

function MOI.get(f::NoObjectiveFilter, attr::MOI.ListOfModelAttributesSet)
    return filter(MOI.get(f.inner, attr)) do a
        return !(a isa MOI.ObjectiveSense || a isa MOI.ObjectiveFunction)
    end
end
function MOI.get(f::NoObjectiveFilter, attr::MOI.NLPBlock)
    block = MOI.get(f.inner, attr)
    return MOI.NLPBlockData(block.constraint_bounds, block.evaluator, false)
end

struct LinearFilter{M<:MOI.ModelLike} <: AbstractModelFilter
    inner::M
end

function MOI.get(f::LinearFilter, attr::MOI.ListOfModelAttributesSet)
    return filter(MOI.get(f.inner, attr)) do a
        return !(a isa MOI.NLPBlock)
    end
end
function MOI.get(f::LinearFilter, attr::MOI.ListOfConstraintTypesPresent)
    return filter(MOI.get(f.inner, attr)) do (F, _)
        return (
            F <: MOI.VariableIndex ||
            F <: MOI.ScalarAffineFunction ||
            F <: MOI.VectorOfVariables ||
            F <: MOI.VectorAffineFunction
        )
    end
end

struct IntegerRelaxation{M<:MOI.ModelLike} <: AbstractModelFilter
    inner::M
end

function MOI.get(f::IntegerRelaxation, attr::MOI.ListOfConstraintTypesPresent)
    return filter(MOI.get(f.inner, attr)) do FS
        S = FS[2]
        return !(S <: MOI.Integer || S <: MOI.ZeroOne)
    end
end

struct FixVariables{T,M<:MOI.ModelLike} <: AbstractModelFilter
    inner::M
    fixed_values::Dict{MOI.VariableIndex,T}
end

function MOI.get(f::FixVariables, attr::MOI.ListOfConstraintTypesPresent)
    ret = MOI.get(f.inner, attr)
    fix = (MOI.VariableIndex, MOI.EqualTo{Float64})
    if !(fix in ret)
        push!(ret, fix)
    end
    return ret
end

function MOI.get(
    f::FixVariables,
    attr::MOI.ListOfConstraintIndices{MOI.VariableIndex,S},
) where {S}
    list = filter(MOI.get(f.inner, attr)) do ci
        return !haskey(f.fixed_values, MOI.VariableIndex(ci.value))
    end
    if S <: MOI.EqualTo
        fix = [
            MOI.ConstraintIndex{MOI.VariableIndex,S}(vi.value) for
            vi in keys(f.fixed_values)
        ]
        list = [list; fix]
    end
    return list
end
function MOI.get(
    f::FixVariables,
    attr::MOI.ConstraintFunction,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,<:MOI.EqualTo},
)
    return MOI.VariableIndex(ci.value)
end
function MOI.get(
    f::FixVariables,
    attr::MOI.ConstraintSet,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,<:MOI.EqualTo},
)
    vi = MOI.VariableIndex(ci.value)
    if haskey(f.fixed_values, vi)
        return MOI.EqualTo(f.fixed_values[vi])
    else
        return MOI.get(f.inner, attr, ci)
    end
end

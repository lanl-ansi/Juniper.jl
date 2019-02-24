# MathOptInterface results
	
function MOI.get(model::Optimizer, ::MOI.TerminationStatus)
	return model.inner.status
end
	
function MOI.get(model::Optimizer, ::MOI.ObjectiveValue)
	if model.inner.status == MOI.OPTIMIZE_NOT_CALLED
	    @error "optimizer not called"
	end
	return model.inner.objval
end
	
function MOI.get(model::Optimizer, ::MOI.ObjectiveBound)
    if model.inner.status == MOI.OPTIMIZE_NOT_CALLED
	    @error "optimizer not called"
	end
	return model.inner.best_bound
end

function MOI.get(model::Optimizer, ::MOI.RelativeGap)
    if model.inner.status == MOI.OPTIMIZE_NOT_CALLED
	    @error "optimizer not called"
	end
	if isnan(model.inner.objval) || isnan(model.inner.best_bound)
		return NaN
	end
	return abs(model.inner.best_bound-model.inner.objval)/abs(model.inner.objval)
end

function MOI.get(model::Optimizer, ::MOI.SolveTime)
    if model.inner.status == MOI.OPTIMIZE_NOT_CALLED
	    @error "optimizer not called"
	end
	return model.inner.soltime
end

function MOI.get(model::Optimizer, ::MOI.VariablePrimal, vi::MOI.VariableIndex)
    if model.inner.status == MOI.OPTIMIZE_NOT_CALLED
	    @error "optimizer not called"
	end
	check_inbounds(model, vi)
	return model.inner.solution[vi.value]
end
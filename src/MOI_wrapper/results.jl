# MathOptInterface results
	
function MOI.get(model::Optimizer, ::MOI.TerminationStatus)
	return model.inner.status
end
	
function MOI.get(model::Optimizer, ::MOI.ObjectiveValue)
    if model.inner.status == MOI.OPTIMIZE_NOT_CALLED
	    @error "optimizer not called"
	end
	return 0.0
end
	
function MOI.get(model::Optimizer, ::MOI.VariablePrimal, vi::MOI.VariableIndex)
	if model.inner.status == MOI.OPTIMIZE_NOT_CALLED
	    @error "optimizer not called"
	end
	check_inbounds(model, vi)
	return model.inner.x[vi.value]
end
export SC_Iterations, MultipleStopCondition

"""
All stopping conditions take a condition that extends the StoppingCondition
abstract type.
"""
abstract type StoppingCondition
end

"""
    SC_Iterations(num)

A stop condition that halts after 'num' iterations.
"""
struct SC_Iterations <: StoppingCondition
    value::Int64
end

"""
    get_SC(constraints...)

Helper function to form a list of stopping conditions.
"""
function get_SC(conditions::StoppingCondition...)::MultipleStopCondition
    return MultipleStopCondition(collect(conditions))
end

"""
    MultipleStopCondition(num)

A stop condition that holds several other stop conditions
"""
struct MultipleStopCondition <: StoppingCondition
    conditions::Vector{StoppingCondition}
end

"""
    compare_MultipleStopCondition(a::MultipleStopCondition,
                                  b::MultipleStopCondition)::Bool

Return true if the conditions within a multiple stop condition type are
valued the same, mainly used for testing.

TODO: make equality independent of vector order.
"""
function compare_MultipleStopCondition(a::MultipleStopCondition,
                                       b::MultipleStopCondition)::Bool
    return a.conditions == b.conditions
end

"""
    stopcondition(v::WorkingVars, checks::MultipleStopCondition)

Checks if any conditions with a MultipleStopCondition type are met.
"""
function stopcondition(model::RAMProblem,
                       checks::MultipleStopCondition
                      )#::Bool
    for c in checks.conditions
		stop, status = stopcondition(model, c)
        if stop
            return true, status
        end
    end
    return false, status
end

"""
    stopcondition(v::WorkingVars, iterations_limit::SC_Iterations)

Checks if the number of iterations has exceeded a maximum.
"""
function stopcondition(model::RAMProblem, 
                       iterations_limit::SC_Iterations
                      )#::Tuple(Bool,internal_termination_conditions)
    return model.iterations >= iterations_limit.value, RAM_ITERATION_LIMIT
end


"""
    check_stopcondition!(model::ModelFormulation, conditions::StoppingCondition)::Bool

Returns true if a condition for stopping has been met. Also updates model's termination
status variable with the value returned by the activated stopping condition. 

Note that the termination variables should be mapped to the MathOptInterface termination 
status codes in MOI_wrapper.jl. If your stopping condition fits an MOI condition that has
not been implemented, then please update the wrapper accordingly.
"""
function check_stopcondition!(model::RAMProblem,
                              conditions::StoppingCondition)::Bool
    stop, status = stopcondition(model, conditions)
	stop && set_termination_status!(model, status)

	return stop
end


function set_termination_status!(model::RAMProblem, status::ram_termination_condition)
    return model.termination_condition = status
end

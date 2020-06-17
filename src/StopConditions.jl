export SC_Iterations

"""
All stopping conditions take a condition that extends the StoppingCondition
abstract type.
"""
abstract type StoppingCondition end

"""
    SC_Iterations(num)

A stop condition that halts after 'num' iterations.
"""
struct SC_Iterations <: StoppingCondition
    value::Int64
end

StopConditionStatus(::SC_Iterations) = ITERATION_LIMIT()



"""
    stopcondition(v::WorkingVars, iterations_limit::SC_Iterations)

Checks if the number of iterations has exceeded a maximum.
"""
function stopcondition(model::RAMProblem, iterations_limit::SC_Iterations)::Bool
    return model.iterations >= iterations_limit.value
end


"""
    check_stopcondition!(model::ModelFormulation, conditions::StoppingCondition)::Bool

Returns true if a condition for stopping has been met. Also updates model's termination
status variable with the value returned by the activated stopping condition. 

Note that the termination variables should be mapped to the MathOptInterface termination 
status codes in MOI_wrapper.jl. If your stopping condition fits an MOI condition that has
not been implemented, then please update the wrapper accordingly.
"""
function check_stopcondition(model::RAMProblem,
                             conditions::Vector{S})::Bool where {S<:StoppingCondition}
    for c in conditions
        stopcondition(model, c) && return true
    end
    return false
end

SetTerminationStatus(m::RAMProblem, c::StoppingCondition) = SetTerminationCondition(m, [c])

function SetTerminationStatus(model::RAMProblem, conditions::Vector{S}) where {S<:StoppingCondition}
    for c in conditions
        !stopcondition(model, c) && continue 
        model.status = StopConditionStatus(c)
        break
    end
end

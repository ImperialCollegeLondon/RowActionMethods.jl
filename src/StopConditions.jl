export IterationStop, TimeStop

"""
All stopping conditions take a condition that extends the StoppingCondition
abstract type.
"""
abstract type StoppingCondition end

SetTerminationStatus(m::RAMProblem, c::StoppingCondition) = SetTerminationCondition(m, [c])

function SetTerminationStatus(model::RAMProblem, conditions::Vector{S}) where {S<:StoppingCondition}
    for c in conditions
        !stopcondition(model, c) && continue
        model.status = StopConditionStatus(c)
        break
    end
end

"""
    _check_stopconditions!(model::ModelFormulation, conditions::StoppingCondition)::Bool

Returns true if a condition for stopping has been met. Also updates model's termination
status variable with the value returned by the activated stopping condition.

Note that the termination variables should be mapped to the MathOptInterface termination
status codes in MOI_wrapper.jl. If your stopping condition fits an MOI condition that has
not been implemented, then please update the wrapper accordingly.
"""
function _check_stopconditions(model::RAMProblem, conditions::Vector{S}
                            )::Bool where {S<:StoppingCondition}

    for c in conditions
        stopcondition(model, c) && return true
    end
    return false
end

"""
    IterationStop(num)

A stop condition that halts after 'num' iterations.
"""
struct IterationStop <: StoppingCondition
    value::Int64
end

StopConditionStatus(::IterationStop) = ITERATION_LIMIT()



"""
    stopcondition(v::WorkingVars, iterations_limit::IterationStop)

Checks if the number of iterations has exceeded a maximum.
"""
function stopcondition(model::RAMProblem, iterations_limit::IterationStop)::Bool
    return model.iterations >= iterations_limit.value
end

mutable struct TimeStop <: StoppingCondition
    running::Bool
    start::Int64
    limit::Int64
    function TimeStop(limit::Int64)
        t = new()
        t.running = false
        t.limit = limit
        return t
    end
end

function stopcondition(model::RAMProblem, time_limit::TimeStop)
    #Start condition if first time
    if time_limit.running == false
        time_limit.running = true
        time_limit.start = floor(time())
        return false
    end

    #Otherwise check if elapsed time is greater than limit
    return floor(time()) - time_limit.start > time_limit.limit
end

StopConditionStatus(::TimeStop) = TIME_LIMIT()

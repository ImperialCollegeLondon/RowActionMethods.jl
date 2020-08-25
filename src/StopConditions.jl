export IterationCondition, TimeCondition

"""
All stopping conditions take a condition that extends the StoppingCondition
abstract type.
"""
abstract type StoppingCondition end

status_code(::StoppingCondition) = OTHER_TERMINATION_CONDITION
status_string(::StoppingCondition) = "Unknown termination criteria reached"
init_stopcondition(::RAMProblem, ::StoppingCondition) = Nothing
test_stopcondition(::RAMProblem, ::StoppingCondition) = error("No stopping condition test function defined")


"""
    _check_stopconditions(model::ModelFormulation, conditions::StoppingCondition)::Bool

Returns true if a condition for stopping has been met. Also updates model's termination
status variable with the value returned by the activated stopping condition.

Note that the termination variables should be mapped to the MathOptInterface termination
status codes in MOI_wrapper.jl. If your stopping condition fits an MOI condition that has
not been implemented, then please update the wrapper accordingly.
"""
function _check_stopconditions(model::RAMProblem, conditions::Vector{S}
                               )::Bool where {S<:StoppingCondition}

    for c in conditions
        if test_stopcondition(model, c)
            model.status = status_code(c)
            model.status_string = status_string(c)
            return false
        end
    end

    return true
end


function _init_stopconditions(model::RAMProblem, conditions::Vector{StoppingCondition})
    for c in conditions
        init_stopcondition(model, c)
    end
end


"""
    IterationCondition(num)

Stop the execution after ``num`` iterations.
"""
struct IterationCondition <: StoppingCondition
    value::Int64

    function IterationCondition(num::Int64)
        if num < 1
            throw( DomainError(num, "Must specify more than one iteration") )
        end

        return new(num)
    end
end

status_code(::IterationCondition) = ITERATION_LIMIT_REACHED
status_string(sc::IterationCondition) = "Reached iteration limit of $(sc.value) iterations"


function test_stopcondition(model::RAMProblem, iterations_limit::IterationCondition)::Bool
    return model.iterations >= iterations_limit.value
end


"""
    TimeCondition(limit::Float64)

Stop the execution of the algorithm after `limit` seconds has elapsed.
"""
mutable struct TimeCondition <: StoppingCondition
    start::Float64
    limit::Float64

    function TimeCondition(limit::Float64)
        if limit <= 0.0
            throw( DomainError(limit, "Time must be postive") )
        end

        return new(0.0, limit)
    end
end


function test_stopcondition(::RAMProblem, time_limit::TimeCondition)
    return ( time() - time_limit.start ) > time_limit.limit
end

function init_stopcondition(::RAMProblem, tc::TimeCondition)
    tc.start = time()
end

status_code(::TimeCondition) = TIME_LIMIT_REACHED
status_string(sc::TimeCondition) = "Reached time limit of $(sc.limit) seconds"

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
    value ::Int64
end

struct MultipleStopCondition <: StoppingCondition
    conditions ::Vector{StoppingCondition}
end

"""
    stopcondition(v::WorkingVars, checks::MultipleStopCondition)

Checks if any conditions with a MultipleStopCondition type are met.
"""
function stopcondition(model::ModelFormulation,
                       checks::MultipleStopCondition
                      )::Bool
    for c in checks.conditions
        if stopcondition(model, c)
            return true
        end
    end
    return false
end

"""
    stopcondition(v::WorkingVars, iterations_limit::SC_Iterations)

Checks if the number of iterations has exceeded a maximum.
"""
function stopcondition(model::ModelFormulation, 
                       iterations_limit::SC_Iterations
                      )::Bool
    return get_iterations(model) >= iterations_limit.value
end


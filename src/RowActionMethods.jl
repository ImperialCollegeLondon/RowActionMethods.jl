module RowActionMethods

#TODO Standardise notation between methods
#=
Author: Edward Stables 
Date: 26-11-2019

For solving QP problems of the form:
Minimize: J = 1/2(x'Ex) + F'x
Subject to the constraints: Mx <= γ
=#

using LinearAlgebra

export iterate_model!, GetModel, buildmodel!, answer, get_SC

include("./Types.jl")
include("./StopConditions.jl")
include("./Hildreth.jl")
include("./ExtendedHildreth.jl")

"""
    get_SC(constraints...)

Helper function to form a list of stopping conditions.
"""
function get_SC(conditions::StoppingCondition...)::MultipleStopCondition
    return MultipleStopCondition(collect(conditions))
end


"""
    iterate_model!(model::ModelFormulation)

Calls iterate_model!(model, condition) with a limit of 32 iterations.
"""
function iterate_model!(model::ModelFormulation)
    conditions = get_SC(SC_Iterations(32)) 
    iterate_model!(model, conditions)
end

"""
    iterate_model!(model::ModelFormulation, conditions::StoppingCondition)

Repeatedly calls model's `iterate!` function until one of the cases
specified in conditions is met. After the condition is met, the final result
will be calculated and returned.
"""
function iterate_model!(model::ModelFormulation, 
                        conditions::StoppingCondition
                       )
    #Returns if unconstrained optimum is valid
    if valid_unconstrained(model)
        set_unconstrained!(model)
        return
    end
    
    #Returns if stopping conditions are already met
    if stopcondition(model, conditions) return end

    #Run iterations until stop conditions are met
    while !stopcondition(model, conditions) 
        iterate!(model)
    end

    #Calculate solution
    resolver!(model)
end

end # module

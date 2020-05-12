mutable struct RAM_Components{T}
    variable_count::Int
    #Maps constraint index var to constraint matrix/value
    constraints::OrderedDict{Int,ConstraintEntry{T}}
    #Tracks largest constraint to ensure a unique new index
    max_constraint_index::Int
    #Track number of constraints
    constraint_count::Int
    termination_condition::ram_termination_condition
    iterations::Int

    function RAM_Components{T}() where {T}
        status = new()
        status.constraints = OrderedDict{Int,ConstraintEntry{T}}()
        status.variable_count = 0
        status.max_constraint_index = 0
        status.constraint_count = 0
        status.termination_condition = RAM_OPTIMIZE_NOT_CALLED
        status.iterations = 0
        return status
    end
end

#TODO this should have a return, and it seems like its use isn't consistent in MOI
function GetModel(model::String)::ModelFormulation
    !haskey(method_mapping, model) && error("Invalid Solver")
    GetModel(method_mapping[model])
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
    #if valid_unconstrained(model)
    #    println("Valid uncons")
    #    set_unconstrained!(model)
    #    return
    #end
    
    #Returns if stopping conditions are already met
    #if stopcondition(model, conditions) return end

    #Run iterations until stop conditions are met
    while !check_stopcondition!(model, conditions) 
        iterate!(model)
        model.status.iterations += 1
    end

    #Calculate solution
    resolver!(model)
end



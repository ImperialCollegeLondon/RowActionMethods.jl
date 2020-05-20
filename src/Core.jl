export GetModel, Setup, Build, Optimize

#TODO this should have a return, and it seems like its use isn't consistent in MOI
function GetModel(model::String)::RAMProblem
    !haskey(method_mapping, model) && 
            throw(ArgumentError("Invalid row action method specified, valid" *
                                " methods are: $(keys(RAM.method_mapping))"))
    return RAMProblem{Float64}(model)
end

ObjectiveType(model::RAMProblem) = ObjectiveType(model.method)

#TODO generalise to other objectives (input vars, and option in else/if)
function Setup(model::RAMProblem{T}, Q::AbstractMatrix, F::AbstractVector) where T
    if ObjectiveType(model.method) == Quadratic()
        model.objective = SparseQuadraticObjective{T}(Q,F)
        model.variable_count = length(F)
    else 
        error("Unsupported objective type")
    end
end

Setup(method::ModelFormulation, model::RAMProblem, Q, F) = Setup(method, Q, F)

Build(model::RAMProblem) = Build(model, model.method)

GetObjective(model::RAMProblem) = GetObjective(model.objective)
GetObjective(obj::SparseQuadraticObjective) = (obj.Q, obj.F)

function Iterate(model::RAMProblem)
    args = [getproperty(model.method, sym) for sym in iterate_args(model.method)]
    if length(args) > 0
        Iterate(model, model.method, args...)
    else
        Iterate(model, model.method)
    end
end

function Resolve(model::RAMProblem)
    args = [getproperty(model.method, sym) for sym in resolve_args(model.method)]
    if length(args) > 0
        Resolve(model, model.method, args...)
    else
        Resolve(model, model.method)
    end
end



"""
    iterate_model!(model::ModelFormulation)

Calls iterate_model!(model, condition) with a limit of 32 iterations.
"""
function Optimize(model::RAMProblem)
    conditions = get_SC(SC_Iterations(32)) 
    Optimize(model, conditions)
end

Optimize(model::RAMProblem, ::Nothing) = Optimize(model)

GetVariables(model::RAMProblem) = GetVariables(model, model.method)

"""
    iterate_model!(model::ModelFormulation, conditions::StoppingCondition)

Repeatedly calls model's `iterate!` function until one of the cases
specified in conditions is met. After the condition is met, the final result
will be calculated and returned.
"""
function Optimize(model::RAMProblem, 
                        conditions::StoppingCondition
                       )
    #Run iterations until stop conditions are met
    while !check_stopcondition!(model, conditions) 
        Iterate(model)
        model.iterations += 1
    end


    #Calculate solution
   Resolve(model)
end


is_empty(model::RAMProblem) = is_empty(model, model.method) && is_model_empty(model)
get_termination_status(model::RAMProblem) = model.termination_condition

ObjectiveValue(model::RAMProblem) = 
    ObjectiveValue(model, ObjectiveType(model.method))

function ObjectiveValue(model::RAMProblem, ::Quadratic)
    B, d = GetObjective(model)
    x = GetVariables(model)

    return 0.5*x'*B*x + x'*d
end


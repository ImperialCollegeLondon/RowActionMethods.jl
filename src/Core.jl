export GetModel, Setup, Optimize, SetThreads, GetVariables, GetObjectiveValue

#TODO this should have a return, and it seems like its use isn't consistent in MOI
function GetModel(model::String; kwargs...)::RAMProblem
    !haskey(method_mapping, model) && 
            throw(ArgumentError("Invalid row action method specified, valid" *
                                " methods are: $(keys(RAM.method_mapping))"))
    return RAMProblem(model; kwargs...)
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

function SetThreads(model::RAMProblem; threads::Bool=true)
    model.threads = threads
end

Setup(method::ModelFormulation, model::RAMProblem, Q, F) = Setup(method, Q, F)

function RunBuild(model::RAMProblem) 
    t1 = time()
    Build(model, model.method)
    model.statistics.BuildTime = time() - t1
end

GetObjective(model::RAMProblem) = GetObjective(model.objective)
GetObjective(obj::SparseQuadraticObjective) = (obj.Q, obj.F)

GetObjectiveFactorised(model::RAMProblem) = GetObjectiveFactorised(model.objective)
GetObjectiveFactorised(obj::SparseQuadraticObjective) = (obj.Qf, obj.F)

Iterate(model::RAMProblem) = Iterate(model.method)

Resolve(model::RAMProblem) = Resolve(model, model.method)

function GetVariables(model::RAMProblem)
    if model.result == nothing
        model.result = GetVariables(model, model.method)
    end
    return model.result
end

"""
    iterate_model!(model::ModelFormulation)

Calls iterate_model!(model, condition) with a limit of 32 iterations.
"""
Optimize(model::RAMProblem) = Optimize(model, [IterationStop(32)])

Optimize(model::RAMProblem, s::StoppingCondition) = Optimize(model, [s])

Optimize(model::RAMProblem, ::Nothing) = Optimize(model)


"""
    Optimize(model::ModelFormulation, conditions::Vector{StoppingCondition})

Repeatedly calls model's `iterate!` function until one of the cases
specified in conditions is met. After the condition is met, the final result
will be calculated and returned.
"""
function Optimize(model::RAMProblem{T}, conditions::Vector{S}) where {T,S<:StoppingCondition}
    
    RunBuild(model)
   
    #Run iterations until stop conditions are met
    #TODO put in checks that the target algorithm supports threading

    t1 = time()
    if !model.threads
        while !check_stopcondition(model, conditions) 
            Iterate(model)
            model.iterations += 1
        end
    else
        #set initial value to the same as that defined in the model
        thread_var = convert(Vector{T}, GetTempVar(model))
        while !check_stopcondition(model, conditions)
            Threads.@threads for i in 1:length(thread_var)
                thread_var[i] = IterateRow(model, i, thread_var) 
            end
            model.iterations += 1
            VarUpdate(model, thread_var)
        end
    end

    SetTerminationStatus(model, conditions)

    #Calculate solution
    #TODO Maybe make this optional, in case the problem will be solved in
    #increments and this doesn't always need to be done?
    Resolve(model)
    model.statistics.OptimizeTime = time() - t1;
end


is_empty(model::RAMProblem) = is_empty(model, model.method) && is_model_empty(model)
get_model_status(model::RAMProblem) = model.status

GetObjectiveValue(model::RAMProblem) = 
    GetObjectiveValue(model, ObjectiveType(model.method))

SupportsDeleteConstraint(model::RAMProblem) = SupportsDeleteConstraint(model.method)
SupportsDeleteConstraint(::ModelFormulation) = false

IterateRow(m::RAMProblem, i::Int, var::Vector{T}) where T = IterateRow(m.method, i, var)
GetTempVar(m::RAMProblem) = GetTempVar(m.method)
VarUpdate(m::RAMProblem{T}, var::Vector{T}) where T = VarUpdate(m.method, var)

function GetObjectiveValue(model::RAMProblem, ::Quadratic)
    B, d = GetObjective(model)
    x = GetVariables(model)

    return 0.5*x'*B*x + x'*d
end


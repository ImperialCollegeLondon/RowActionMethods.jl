export GetModel, Optimize, SetThreads, GetVariables, GetObjectiveValue

"""
    GetModel(model::String; kwargs...)::RAMProblem

Return problem definition configured with the specified method.

kwargs are passed directly to the method's constructor, see the method for 
documentation on available options.

This function queries `RAM.method_mapping` to find valid methods.
"""
function GetModel(method::String; kwargs...)::RAMProblem
    !haskey(method_mapping, method) && 
            throw(ArgumentError("Invalid row action method specified, valid" *
                                " methods are: $(keys(RAM.method_mapping))"))
    return RAMProblem(method; kwargs...)
end

ObjectiveType(model::RAMProblem) = ObjectiveType(model.method)

"""
    SetThreads(model::RAMProblem; threads::Bool=true)

Enable/disable threading for `model`.

This only has an effect if the method used by `model` has implemented 
a threaded variation of its algorithm.

Julia will set the number of threads to that defined by the
`JULIA_NUM_THREADS` environment variable.
"""
function SetThreads(model::RAMProblem; threads::Bool=true)
    model.threads = threads
end

"""
    GetObjective(model::RAMProblem)::Tuple

Return the components of the objective function set for `model`. This function
is implemented for each of the objective types.
"""
GetObjective(model::RAMProblem)::Tuple = GetObjective(model.objective)

"""
    GetObjective(obj::SparseQuadraticObjective)::Tuple

Return the `Q` matrix and `F` vector of the objective.
"""
GetObjective(obj::SparseQuadraticObjective)::Tuple = (obj.Q, obj.F)

"""
    GetObjectiveFactorised(model::RAMProblem)::Tuple

As `GetObjective` but for factorised matrices.
"""
GetObjectiveFactorised(model::RAMProblem)::Tuple = GetObjectiveFactorised(model.objective)

"""
    GetObjectiveFactorised(obj::SparseQuadraticObjective)::Tuple

Return the factorised `Q` matrix (Cholseky)and `F` vector of the objective. Note that
the factorised `Q` matrix is not sparse.
"""
GetObjectiveFactorised(obj::SparseQuadraticObjective)::Tuple = (obj.Qf, obj.F)

"""
    is_empty(model::RAMProblem)::Bool

Return true if a model has not had any parameters set.
"""
is_empty(model::RAMProblem)::Bool = is_empty(model, model.method) && is_model_empty(model)

"""
    get_model_status(model::RAMProblem)

Getter function for `model.status`.
"""
get_model_status(model::RAMProblem) = model.status

#Modification defaults
SupportsDeleteConstraint(model::RAMProblem) = SupportsDeleteConstraint(model.method)
SupportsDeleteConstraint(::ModelFormulation) = false

#Threaded function defaults
IterateRow(m::RAMProblem, i::Int, var::Vector{T}) where T = IterateRow(m.method, i, var)
GetTempVar(m::RAMProblem) = GetTempVar(m.method)
VarUpdate(m::RAMProblem{T}, var::Vector{T}) where T = VarUpdate(m.method, var)

"""
    GetObjectiveValue(model::RAMProblem{T,F})::T
    
Return the evaluated objective for the current solution.
"""
(GetObjectiveValue(model::RAMProblem{T,F})::T) where {T,F}=
    GetObjectiveValue(model, ObjectiveType(model.method))


"""
    GetObjectiveValue(model::RAMProblem, ::Quadratic)

Evaluate the cost of the objective function for a quadratic objective.
"""
function GetObjectiveValue(model::RAMProblem, ::Quadratic)
    B, d = GetObjective(model)
    x = GetVariables(model)

    return 0.5*x'*B*x + x'*d
end

SetObjective(model::RAMProblem, args...) = 
    SetObjective(model, ObjectiveType(model.method), args...)

"""
    SetObjective(model::RAMProblem{T}, ::Quadratic, Q::AbstractMatrix, F::AbstractVector) where T

Set the objective of `model` as a quadratic type.
"""
function SetObjective(model::RAMProblem{T,F}, ::Quadratic, Q::AbstractMatrix, V::AbstractVector) where {T,F}
    model.objective = SparseQuadraticObjective{T}(Q,V)
    model.variable_count = length(V)
end

#Build wrapper to set timings
function RunBuild(model::RAMProblem) 
    t1 = time()
    Build(model, model.method)
    model.statistics.BuildTime = time() - t1
end

Optimize(model::RAMProblem, ::Nothing) = Optimize(model)

"""
    Optimize(model::RAMProblem)

Iterate `model.method` algorithm until a stopping condition has been met, then
set the termination status and calculate the primal solution.

If configured as single threaded, `Iterate` is called on `model`, with each call
to `Iterate` counting as a single iteration.

If configured as multi threaded then `IterateRow` is called for each index in the
variable returned by `GetTempVar`. Calls to `IterateRow` are distributed amongst all
available threads.

Configures the algorithm to terminate after 32 iterations.
"""
Optimize(model::RAMProblem) = Optimize(model, [IterationStop(32)])

"""
    Optimize(model::RAMProblem, s::StoppingCondition)

As [`Optimize`](@ref Optimize(::RAM.RAMProblem)) but takes a single stopping condition to end on rather than
the default.
"""
Optimize(model::RAMProblem, s::StoppingCondition) = Optimize(model, [s])

"""
    Optimize(model::RAMProblem{T,F}, conditions::Vector{S}) where {T,F,S<:StoppingCondition}

As [`Optimize`](@ref Optimize(::RAM.RAMProblem)) but takes a vector of single stopping conditions to end on. 
"""
function Optimize(model::RAMProblem{T,F}, conditions::Vector{S}) where {T,F,S<:StoppingCondition}
    
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

#Wrapper
Iterate(model::RAMProblem) = Iterate(model.method)

#Wrapper
Resolve(model::RAMProblem) = Resolve(model, model.method)

#Wrapper
function GetVariables(model::RAMProblem)
    if model.result == nothing
        model.result = GetVariables(model, model.method)
    end
    return model.result
end



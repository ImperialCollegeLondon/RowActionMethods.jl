export GetModel

mutable struct RAMProblem{T}
    #== Variables ==#
    variable_count::Int
    
    #== Constraints ==#
    #Maps constraint index var to constraint matrix/value
    constraints::OrderedDict{Int,ConstraintEntry{T}}
    #Tracks largest constraint to ensure a unique new index
    max_constraint_index::Int
    #Track number of constraints
    constraint_count::Int

    #== Problem Description ==#
    #TODO Cache most recent indexes possibly?
    SparseMatrices::Vector{SparseMatrixCSC{T,Int64}}
    SparseVectors::Vector{SparseVector{T,Int64}}

    #== General ==#
    termination_condition::ram_termination_condition
    iterations::Int
    method::ModelFormulation

    function RAMProblem{T}(model::String) where T
        p = new()

        #== Constraints ==#
        p.constraints = OrderedDict{Int,ConstraintEntry{T}}()
        p.max_constraint_index = 0
        p.constraint_count = 0

        #== Problem ==#
        p.SparseMatrices = []
        p.SparseVectors = []

        p.variable_count = 0
        p.termination_condition = RAM_OPTIMIZE_NOT_CALLED
        p.iterations = 0

        p.method = method_mapping[model]{T}()
        return p 
    end
end

struct SparseMatrixIndex
    value::Int64
end

struct SparseVectorIndex
    value::Int64
end

function RegisterSparse(model::RAMProblem{T}, data::Matrix)::SparseMatrixIndex where T
    push!(model.SparseMatrices, convert(Matrix{T}, data))
    return SparseMatrixIndex(length(model.SparseMatrices) - 1)
end

function RegisterSparse(model::RAMProblem{T}, data::Vector)::SparseVectorIndex where T
    push!(model.SparseVectors, convert(Vector{T}, data))
    return SparseVectorIndex(length(model.SparseVectors) - 1)
end

GetSparse(model::RAMProblem, index::SparseMatrixIndex)::SparseMatrixCSC = model.SparseMatrices[index.value]

GetSparse(model::RAMProblem, index::SparseVectorIndex)::SparseVector = model.SparseVectors[index.value]

#TODO this should have a return, and it seems like its use isn't consistent in MOI
function GetModel(model::String)::RAMProblem
    !haskey(method_mapping, model) && error("Invalid Solver")
    return RAMProblem{Float64}(model)
end

#TODO generalise to other objectives (input vars, and option in else/if)
function Setup(model::RAMProblem, Q::Matrix, F::Vector, variables::Int)
    if ObjectiveType(model.method) == Quadratic()
        Q = RegisterSparse(model, Q)  
        F = RegisterSparse(model, F)
        model.variable_count = variables
        Setup(model.method, model, Q, F)
    else 
        error("Unsupported objective type")
    end
end

Setup(method::ModelFormulation, model::RAMProblem, Q, F) = Setup(method, Q, F)

Build(model::RAMProblem) = Build(model.method, model)

function Iterate(model::RAMProblem)
    args = [getproperty(model.method, sym) for sym in iterate_args(model.method)]
    Iterate(model.method, args...)
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
    #Run iterations until stop conditions are met
    while !check_stopcondition!(model, conditions) 
        Iterate(model)
        model.status.iterations += 1
    end

    #Calculate solution
    resolver!(model)
end



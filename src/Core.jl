export GetModel


struct SparseMatrixIndex
    value::Int64
end

struct SparseVectorIndex
    value::Int64
end

function RegisterSparse(model::RAMProblem{T}, data::Matrix)::SparseMatrixIndex where T
    push!(model.SparseMatrices, convert(Matrix{T}, data))
    return SparseMatrixIndex(length(model.SparseMatrices))
end

function RegisterSparse(model::RAMProblem{T}, data::Vector)::SparseVectorIndex where T
    push!(model.SparseVectors, convert(Vector{T}, data))
    return SparseVectorIndex(length(model.SparseVectors))
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
    args = [GetSparse(model, getproperty(model.method, sym)) for sym in iterate_args(model.method)]
    Iterate(model, model.method, args...)
end

function Resolve(model::RAMProblem)
    args = [GetSparse(model, getproperty(model.method, sym)) for sym in resolve_args(model.method)]
    Resolve(model, model.method, args...)
end

"""
    iterate_model!(model::ModelFormulation)

Calls iterate_model!(model, condition) with a limit of 32 iterations.
"""
function Optimize(model::RAMProblem)
    conditions = get_SC(SC_Iterations(32)) 
    Optimize(model, conditions)
end

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



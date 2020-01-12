import MathOptInterface

const MOI = MathOptInterface
const RAM = RowActionMethods

mutable struct Optimizer <: MOI.AbstractOptimizer

    method::RAM.RowActionMethod
    inner_model::RAM.ModelFormulation
    sense::MOI.OptimizationSense
    Name::String

    variable_count::Int

    function Optimizer(method::String;kwargs...)
        model = new()

        if !haskey(RAM.method_mapping, method) 
            throw(ArgumentError("Invalid row action method specified, valid" *
                                " methods are: $(keys(RAM.method_mapping))"))
        end
        
        model.Name = "RowActionMethods-$method"
        model.method = RAM.method_mapping[method]
        model.inner_model = GetModel(model.method)

        for (key, val) in kwargs
            MOI.set(model, MOI.RawParameter(String(key)), val)
        end

        MOI.empty!(model)

        return model
    end
end

function MOI.supports(::Optimizer, 
                      ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}})
    return true
end

function MOI.empty!(model)
    model.variable_count = 0
end

function MOI.add_variables(model::Optimizer, n::Int)::Vector{MOI.VariableIndex}
    model.variable_count = n
    
    return [MOI.VariableIndex(i) for i in 1:n]
end

function MOI.set(model::Optimizer, 
                 ::MOI.ObjectiveFunction{F}, 
                 val::F) where {F <: MOI.ScalarQuadraticFunction{Float64}} 

    Q = zeros(model.variable_count, model.variable_count)
    a = zeros(model.variable_count)
    
    for t in val.quadratic_terms
        Q[t.variable_index_1.value, t.variable_index_2.value] = t.coefficient
    end

    for t in val.affine_terms
        a[t.variable_index.value] = t.coefficient
    end

    setfunction!(model.inner_model, Q, a)
end

#=
function MOI.add_constraint(model::Optimizer,
                            func_list::Array{MOI.ScalarAffineFunction{T}, 1},
                            lim_list::Array{S{T}, 1}) where {T, S<:MOI.AbstractScalarSet}
    if sizeof(func_list)[1] != sizeof(lim_list[1])
        throw(BoundsError("The number of constaint equations must match the number of constraint inequalities."))
    end

    for constraint in zip(func_list, lim_list)
        MOI.add_constraint(model, constraint[1], constraint[2])
    end
end
=#

function MOI.add_constraint(model::Optimizer, 
                            func::MOI.ScalarAffineFunction, 
                            lim::MOI.LessThan)
    #TODO: Raise error on non-zero valued constant in func
    constraint_function = zeros(model.variable_count)
    for t in func.terms
        constraint_function[t.variable_index.value] = t.coefficient
    end

    setconstraints!(model.inner_model, constraint_function, lim.upper)
end

function MOI.add_constraint(model::Optimizer,
                            func::MOI.ScalarAffineFunction,
                            lim::MOI.GreaterThan)
    MOI.add_constraint(model, func, MOI.LessThan(-lim))
end

function MOI.is_empty(model::Optimizer)::Bool
    return RAM.is_empty(model.inner_model)
end

function MOI.optimize!(model::Optimizer)
    RAM.buildmodel!(model.inner_model)
    RAM.iterate_model!(model.inner_model)
end

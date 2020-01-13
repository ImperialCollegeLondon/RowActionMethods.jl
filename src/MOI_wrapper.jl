import MathOptInterface

const MOI = MathOptInterface
const MOIU = MOI.Utilities
const RAM = RowActionMethods

mutable struct Optimizer <: MOI.AbstractOptimizer

    method::RAM.RowActionMethod
    inner_model::RAM.ModelFormulation
    sense::MOI.OptimizationSense
    name::String

    variable_count::Int

    silent::Bool #True means nothing should be printed

    function Optimizer(method::String;kwargs...)
        model = new()

        if !haskey(RAM.method_mapping, method) 
            throw(ArgumentError("Invalid row action method specified, valid" *
                                " methods are: $(keys(RAM.method_mapping))"))
        end
        
        model.name = "RowActionMethods-$method"
        model.method = RAM.method_mapping[method]
        model.inner_model = GetModel(model.method)
        model.silent= false

        for (key, val) in kwargs
            MOI.set(model, MOI.RawParameter(String(key)), val)
        end

        MOI.empty!(model)

        return model
    end
end

#= Model Actions =#
function MOI.optimize!(model::Optimizer)
    RAM.buildmodel!(model.inner_model)
    RAM.iterate_model!(model.inner_model)
end

#= Model Status =#
function MOI.supports(::Optimizer, 
                      ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}})
    return true
end

function MOI.is_empty(model::Optimizer)::Bool
    return RAM.is_empty(model.inner_model)
end

#= Model special set functions =#
function MOI.empty!(model)
    model.variable_count = 0
end

#= Variables =#
function MOI.add_variables(model::Optimizer, n::Int)::Vector{MOI.VariableIndex}
    model.variable_count = n
    
    return [MOI.VariableIndex(i) for i in 1:n]
end

#= Constraints =#
function MOI.add_constraint(model::Optimizer, 
                            func::MOI.ScalarAffineFunction{T}, 
                            lim::MOI.LessThan{T}) where T
    #TODO: Raise error on non-zero valued constant in func
    constraint_function = zeros(model.variable_count)
    for t in func.terms
        constraint_function[t.variable_index.value] = t.coefficient
    end

    RAM.setconstraints!(model.inner_model, constraint_function, lim.upper)
end

function MOI.add_constraint(model::Optimizer,
                            func::MOI.ScalarAffineFunction,
                            lim::MOI.GreaterThan)
    MOI.add_constraint(model, func, MOI.LessThan(-lim.lower))
end

#= MOI.set functions =#
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

function MOI.set(model::Optimizer, ::MOI.Silent, val::Bool)
    model.silent = val
end

#= MOI.get functions =#
function MOI.get(model::Optimizer, ::MOI.SolverName)::String
    return model.name
end

#= MOI Copy Functions =#
function MOI.copy_to(dest::Optimizer, src::MOI.ModelLike; kws...)
    return MOI.Utilities.automatic_copy_to(dest, src; kws...)
end

MOIU.supports_default_copy_to(model::Optimizer, copy_names::Bool) = true



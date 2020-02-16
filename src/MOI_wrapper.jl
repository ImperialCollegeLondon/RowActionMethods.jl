using DataStructures
import MathOptInterface

const MOI = MathOptInterface
const MOIU = MOI.Utilities
const MOIB = MOI.Bridges
const RAM = RowActionMethods

mutable struct Optimizer <: MOI.AbstractOptimizer

    method::RAM.RowActionMethod
    inner_model::RAM.ModelFormulation
    sense::MOI.OptimizationSense
    name::String

    variable_count::Int
    constraints::Vector{MOI.ConstraintIndex}

    silent::Bool #True means nothing should be printed

    function Optimizer(method::String;kwargs...)
        model = new()

        if !haskey(RAM.method_mapping, method) 
            throw(ArgumentError("Invalid row action method specified, valid" *
                                " methods are: $(keys(RAM.method_mapping))"))
        end
        
        model.name = "RowActionMethods-$method"
        model.method = RAM.method_mapping[method]
        model.silent = false
        model.variable_count = 0
        model.constraints = Vector{MOI.ConstraintIndex}()
                                  

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
                      ::MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}})
    return true
end

#Note that this is a temporary measure, there is no reason to not support max but
#at the time of writing it requires refactoring code that will break even more tests
MOI.supports(::Optimizer, ::MOI.ObjectiveSense) = true

function MOI.is_empty(model::Optimizer)::Bool
    return RAM.is_empty(model.inner_model)
end

#= Model special set functions =#
function MOI.empty!(model)
    model.inner_model = GetModel(model.method)
    model.variable_count = 0
    model.constraints = []
end

#= Variables =#
function MOI.add_variable(model::Optimizer)::MOI.VariableIndex
    model.variable_count += 1
    return MOI.VariableIndex(model.variable_count)
end

function MOI.add_variables(model::Optimizer, n::Int)::Vector{MOI.VariableIndex}
    return [MOI.add_variable(model) for i in 1:n]
end

function MOI.delete(model::Optimizer, index::MOI.VariableIndex)
    #TODO add support for deletion when constraints exist
    #TODO understand what affect this has on objective
    if length(model.constraints) > 0
        error("Cannot delete variables when constraints exist.")
    end
    RAM.delete_variable!(model.inner_model)
end

#= Constraints =#
#TODO: Add ability to add non lessthan constraints, probably needs indicator 
#function in solvers?
function MOI.add_constraint(model::Optimizer, 
                            func::MOI.ScalarAffineFunction{T}, 
                            lim::MOI.LessThan{T})::MOI.ConstraintIndex where T
    #TODO: Raise error on non-zero valued constant in func
    constraint_function = zeros(model.variable_count)
    for t in func.terms
        constraint_function[t.variable_index.value] = t.coefficient
    end

    index = RAM.addconstraint!(model.inner_model, constraint_function, lim.upper)
    constraint_index = MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64}, 
                                           MOI.LessThan{Float64}}(index) 
    push!(model.constraints, constraint_index)
    return constraint_index
end

function MOI.supports_constraint(model::Optimizer,
                                 ::Type{MOI.ScalarAffineFunction{T}},
                                 ::Type{MOI.LessThan{T}})::Bool where T
    return true
end

function MOI.supports_constraint(model::Optimizer,
                                 ::Type{MOI.VectorOfVariables},
                                 ::Type{MOI.Reals})::Bool 
    return true
end

function MOI.delete(model::Optimizer, index::MOI.ConstraintIndex{F,S}) where {F,S}
    RAM.delete_constraint!(model.inner_model, index.value)
    #TODO alternative way to remove entry efficiently?
    setdiff!(model.constraints, [index])
end

function MOI.modify(model::Optimizer, 
                    con_index::MOI.ConstraintIndex{F,S},
                    change::MOI.ScalarCoefficientChange{Float64}) where {F,S}
    RAM.edit_constraint_coefficient!(model.inner_model, 
                                     con_index.value, 
                                     change.variable.value, 
                                     change.new_coefficient)
end

function MOI.modify(model::Optimizer,
                    con_index::MOI.ConstraintIndex{F,S},
                    change::MOI.ScalarConstantChange{Float64}) where {F,S}
    RAM.edit_constraint_constant!(model.inner_model,
                                  con_index.value, 
                                  change.new_constant)
end

#Currently not supporting in-place modification, rebuilds the model each time
function MOI.set(model::Optimizer, ::MOI.ConstraintSet, c::MOI.ConstraintIndex{F,MOI.LessThan}, set::MOI.LessThan) where {F, S}
    MOI.modify(model, c, MOI.ScalarConstantChange(c.upper))    
end

#= MOI.set functions =#
function MOI.set(model::Optimizer, 
                 ::MOI.ObjectiveFunction{F}, 
                 val::F) where {F <: MOI.ScalarQuadraticFunction{Float64}} 

    num_vars = model.variable_count
    Q = zeros(num_vars, num_vars)
    a = zeros(num_vars)
    
    for t in val.quadratic_terms
        Q[t.variable_index_1.value, t.variable_index_2.value] = t.coefficient
        Q[t.variable_index_2.value, t.variable_index_1.value] = t.coefficient
    end

    for t in val.affine_terms
        a[t.variable_index.value] = t.coefficient
    end

    setobjective!(model.inner_model, Q, a, num_vars)
end

function MOI.set(model::Optimizer, ::MOI.Silent, val::Bool)
    model.silent = val
end

function MOI.set(model::Optimizer, ::MOI.ObjectiveSense, val::MOI.OptimizationSense)
    model.sense = val
end


#= MOI.get functions =#
function MOI.get(model::Optimizer, ::MOI.SolverName)::String
    return model.name
end

function MOI.get(model::Optimizer, ::MOI.NumberOfVariables)::Int
    return model.variable_count
end

function MOI.get(model::Optimizer, 
                 ::MOI.NumberOfConstraints{MOI.ScalarAffineFunction{T},
                                           MOI.LessThan{T}})::Int where T
    #FIXME update with type access not dict access
    return model.constraint["LessThan"]
end

#Maps internal RAM termination conditions to MOI equivalents
RAM_MOI_Termination_Map = Dict(
RAM.RAM_OPTIMIZE_NOT_CALLED => MOI.OPTIMIZE_NOT_CALLED,
RAM.RAM_OPTIMAL             => MOI.OPTIMAL,
RAM.RAM_INFEASIBLE          => MOI.INFEASIBLE,
RAM.RAM_ITERATION_LIMIT     => MOI.ITERATION_LIMIT,  
RAM.RAM_TIME_LIMIT          => MOI.TIME_LIMIT
)

function MOI.get(model::Optimizer, ::MOI.TerminationStatus)
    return RAM_MOI_Termination_Map[RAM.get_termination_status(model.inner_model)]
end

#function MOI.get(model::Optimizer, ::MOI.ObjectiveValue)
#    sense = model.sense == MOI.MAX_SENSE ? -1 : 1
#    return sense * RAM.answer(model.inner_model)
#end

function MOI.get(model::Optimizer, var::MOI.VariablePrimal, vi::MOI.VariableIndex)
    sense = model.sense == MOI.MAX_SENSE ? -1 : 1
    return sense * RAM.answer(model.inner_model)[vi.value]
end

#= MOI Copy Functions =#
function MOI.copy_to(dest::Optimizer, src::MOI.ModelLike; kws...)
    return MOI.Utilities.automatic_copy_to(dest, src; kws...)
end

MOIU.supports_default_copy_to(model::Optimizer, copy_names::Bool) = true


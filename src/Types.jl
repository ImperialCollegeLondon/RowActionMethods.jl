using SuiteSparse
"""
All variables used in the actual computation should be stored here. Use of dictionary gives flexibility to contents, leaves the door open to changing the iterator between iterations.
"""


"""
Type used to find iteration functions and building initial
model. Extend in file for each algorithm type.
"""
abstract type RowActionMethod end

"""
Type used as supertype of specific model formulations used by 
different algorithms.
"""
abstract type ModelFormulation end

abstract type AbstractObjectiveType end
struct Quadratic <: AbstractObjectiveType end 
struct Linear <: AbstractObjectiveType end 
ObjectiveType(m::ModelFormulation) = error("$(typeof(m)) should define an objective function type")

abstract type AbstractObjective end

#TODO make Qf a sparse type
struct SparseQuadraticObjective{T} <: AbstractObjective
    Q::SparseMatrixCSC{T}
    #Qf::SuiteSparse.CHOLMOD.Factor{T}
    Qf::Cholesky{T}
    F::SparseVector{T}
    function SparseQuadraticObjective{T}(Q,F) where T
        Qf = cholesky(Q)
        Q = sparse(Q)
        F = sparse(F)
        return new(Q, Qf, F)
    end
end

#TODO Add sparse vectors in Constraint entry
mutable struct ConstraintEntry{T}
    func::Vector{T}
    lim::T
end

@enum(ram_termination_condition,
RAM_OPTIMIZE_NOT_CALLED,
RAM_OPTIMAL, RAM_INFEASIBLE,
RAM_ITERATION_LIMIT, RAM_TIME_LIMIT)

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
    objective::AbstractObjective

    #== General ==#
    termination_condition::ram_termination_condition
    iterations::Int
    method::ModelFormulation

    function RAMProblem{T}(model::String; kwargs...) where T
        p = new()

        #== Constraints ==#
        p.constraints = OrderedDict{Int,ConstraintEntry{T}}()
        p.max_constraint_index = 0
        p.constraint_count = 0

        p.variable_count = 0
        p.termination_condition = RAM_OPTIMIZE_NOT_CALLED
        p.iterations = 0

        p.method = method_mapping[model]{T}(;kwargs...)
        return p 
    end
end


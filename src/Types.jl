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

struct SparseQuadraticObjective{T} <: AbstractObjective
    Q::SparseMatrixCSC{T}
    #Qf::SuiteSparse.CHOLMOD.Factor{T}
    Qf::Cholesky{T}
    F::SparseVector{T}
    function SparseQuadraticObjective{T}(Q,F) where T
        #TODO this formulation needs to be changed
        #currently Q must be dense, but this is memory
        #inefficient and wasteful if sparse is provided 
        #in the first place
        Qf = cholesky(Q)
        Q = sparse(Q)
        F = sparse(F)
        return new(Q, Qf, F)
    end
end

mutable struct Statistics{T}
    BuildTime::T
    OptimizeTime::T
    Statistics{T}() where T = new(0.0, 0.0)
end


#TODO build out type heirarchy
mutable struct ConstraintEntry{T}
    func::SparseVector{T}
    lim::T
end

abstract type AbstractStatus end

struct OPTIMIZE_NOT_CALLED              <: AbstractStatus end
struct OPTIMAL                          <: AbstractStatus end
struct INFEASIBLE                       <: AbstractStatus end
struct ITERATION_LIMIT                  <: AbstractStatus end
struct TIME_LIMIT                       <: AbstractStatus end
struct UNKNOWN_TERMINATION_CONDITION    <: AbstractStatus end

mutable struct RAMProblem{T,F}
    #== Variables ==#
    variable_count::F
    
    #== Constraints ==#
    #Maps constraint index to actual vector index
    constraint_indexes::Dict{F, F}
    #Maps constraint index to constraint vector
    constraints::OrderedDict{F,ConstraintEntry{T}}
    #Tracks largest constraint to ensure a unique new index
    max_constraint_index::F
    #Track number of constraints
    constraint_count::F

    #== Problem Description ==#
    objective::AbstractObjective
    result::Union{SparseVector{T},Nothing}

    #== General ==#
    status::AbstractStatus
    iterations::F
    method::ModelFormulation

    #== Threading ==#
    threads::Bool
    
    statistics::Statistics{T}

    RAMProblem(model::String; kwargs...) = RAMProblem{Float64, Int64}(model; kwargs...)
    function RAMProblem{T,F}(model::String; kwargs...) where {T,F}
        p = new()

        #== Constraints ==#
        p.constraint_indexes = Dict{F,F}()
        p.constraints = OrderedDict{F,ConstraintEntry{T}}()
        p.max_constraint_index = 0
        p.constraint_count = 0

        p.variable_count = 0
        p.status = OPTIMIZE_NOT_CALLED()
        p.iterations = 0
        p.result = nothing
        p.threads = false

        p.statistics = Statistics{T}()

        p.method = method_mapping[model]{T}(;kwargs...)
        return p 
    end
end


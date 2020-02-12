using DataStructures
import Base.==
export Hildreth

struct Hildreth <: RowActionMethod end

#TODO Implement vectorisation
"""
    HildrethModel(E, F, M, γ, H K, ucSoln, Soln, E_fact, workingvars, options)

A model for solving the problem with Hildreth's original method.
For the problem: 
min ½(x'Ex)+F'x s.t. Mx≦γ

ucSoln - Unconstrained solution of the problem\n
Soln - calculated, constrained solution of the problem\n
E_fact - factorised representation of E\n
workingvars - a Dict of values adjusted between iterations\n
options - [currently unused] specific values used for setting solver-specific values
"""
#TODO change to parametric type
mutable struct HildrethModel <: ModelFormulation
    variable_count::Int
    #QP matrix
    E::Array{Float64}
    #QP vector
    F::Vector{Float64}
    #Maps constraint index var to constraint matrix/value
    constraints::OrderedDict{Int,Pair{Vector{Float64},Float64}}
    #Tracks largest constraint to ensure a unique new index
    max_constraint_index::Int
    #Track number of constraints
    constraint_count::Int 
    #M::Array{Float64}
    #γ::Vector{Float64}
    H::Array{Float64}
    K::Vector{Float64}
    ucSoln::Vector{Float64}
    Soln::Vector{Float64}
    E_fact::Union{Bidiagonal,Factorization,Array,Diagonal}
    workingvars::Dict{String,Any}
    termination_condition::internal_termination_conditions

    function HildrethModel()
        model = new()
        model.workingvars = Dict("iterations" => 0)
        model.termination_condition = RAM_OPTIMIZE_NOT_CALLED

        model.constraints = OrderedDict{Int,Pair{Vector{Float64},Float64}}()

        model.variable_count = 0
        model.constraint_count = 0
        model.max_constraint_index = 0
        return model
    end
end

function get_termination_status(model::HildrethModel)::internal_termination_conditions
    return model.termination_condition
end

function set_termination_status!(model::HildrethModel, status::internal_termination_conditions)
    return model.termination_condition = status
end

"""
    GetModel(::Hildreth)::HildrethModel
    
Returns a seleton model of the problem for solving with Hildreth's orignal
method.
"""
function GetModel(::Hildreth)::HildrethModel
    return HildrethModel()
end

"""
    iterate!(model::HildrethModel)

Performs one iteration of the algorithm. Updates λ as it progresses.

Treats the entire summation as a calculation of H_i * λ, then subtracts the 
contribution of the currently considered λ. 
"""
function iterate!(model::HildrethModel)
    λ =  model.workingvars["λ"]
    model.workingvars["λ_old"] = copy(λ)
    H = model.H
    K = model.K
    
    for (i,l) in enumerate(λ)
        w = (H[i:i,:] * λ)[1] - H[i,i] * l
        w += K[i]
        w /= -H[i,i]
        λ[i] = max(0, w)
    end

    model.workingvars["iterations"] += 1
end

"""
    setobjective!(model::HildrethModel, E::Array{T, 2}, F::Vector{T}) where T

Sets the objective function for the problem. Hildreth's algorithm requires that 
E is positive definite, and the dimensions of both to match the number of variables
defined in the problem (if using the MOI/JuMP interface this should be enforced).

Currently store all matrices in a standard format. As row action matrices are
intended for sparse matrices, this should probably be the type used. In addition,
we have a requirement for positive definite matrices, this is assumed in insertion
and should be checked. 
"""
function setobjective!(model::HildrethModel, E::Array{T, 2}, F::Vector{T}, num_vars) where T
    model.E = E
    model.F = F
    model.variable_count = num_vars
end

"""
    setconstriant!(model::HildrethModel, M_row::Vector{T}, lim::T) where T

Adds a less-than constraint to the model.

Returns a UID value that the solver can use to map to the value if modifying,
viewing, or deleting the values. UID is based off the number of constraints
that have ever been added, not the current number.
"""
function addconstraint!(model::HildrethModel, M_row::Vector{T}, lim::T)::Int where T
    !validconstraint(model, M_row, lim) && error("Invalid constraint")
    
    #Ensures unique constraint index
    new_index = model.max_constraint_index + 1

    push!(model.constraints, new_index => M_row => lim)

    #Update largest index
    model.max_constraint_index = new_index
    model.constraint_count += 1
    return new_index
end

"""
    validconstraint(model::HildrethModel, row::Vector{T}, lim::T)::Bool where T

Return a bool to indicate if a constraint is currently valid for the model.

Checks if the number of entries is equal to the number of registered variables.
"""
function validconstraint(model::HildrethModel, row::Vector{T}, lim::T)::Bool where T
    return size(row)[1] == model.variable_count
end

"""
    extendcontraints(model::HildrethModel)
    
Append a new value to the end of each constraint in model. Assumes that a new
variable being added is added at the end.
"""
function extendcontraints(model::HildrethModel)
    #TODO add type parameters for appending
    for (_,(r,_)) in model.constraints
        append!(r, 0.0)
    end
end

"""
    shrinkconstraints(model::HildrethModel, index::Int)
    
Removes the entry in each constraint at index.
"""
function shrinkconstraints(model::HildrethModel, index::Int)
    for (_,(r,_)) in model.constraints
        deleteat!(r, index)
    end
end

"""
    get_constraintmatrix(model::HildrethModel)::Matrix(Float64,2)

Forms the stored constraints into a transposed matrix form,
ie the value Mᵀ is returned.
"""
#TODO: Add type parameterisation
function get_constraintmatrix(model::HildrethModel)::Array{Float64,2}
    a = hcat([i for (_,(i,_)) in model.constraints]...)
    println(a)
    return a
end

#TODO: Add type parameterisation
function get_constraintvector(model::HildrethModel)::Vector{Float64}
    return [j for (_,j) in values(model.constraints)]
end

"""
    buildmodel(model::HildrethModel)

Builds the internal variables based on problem specification
"""
function buildmodel!(model::HildrethModel)
    Mt = get_constraintmatrix(model)
    M = Mt'
    γ = get_constraintvector(model)

    model.E_fact = cholesky(model.E)
    model.H = M * (model.E_fact\Mt)
    model.K = γ + (M * (model.E_fact\model.F))
    model.ucSoln = -(model.E_fact\model.F)
    model.workingvars["λ"] = zeros(model.constraint_count)
    model.workingvars["λ_old"] = zeros(model.constraint_count)
end

"""
    ==(a::HildrethModel, b::HildrethModel)::Bool

Checks equality of two hildreth model structs. Does not check for object
equality, only value equality.
"""
function ==(a::HildrethModel, b::HildrethModel)::Bool
    return a.E == b.E &&
           a.F == b.F &&
           a.M == b.M &&
           a.γ == b.γ &&
           a.E_fact == b.E_fact &&
           a.H == b.H &&
           a.K == b.K &&
           a.ucSoln == b.ucSoln &&
           a.workingvars == b.workingvars
end

"""
    valid_unconstrained(model::HildrethModel)::Bool

Returns true if the initial minimum is within the constraints.
"""
function valid_unconstrained(model::HildrethModel)::Bool
    valid = model.M * model.ucSoln - model.γ
    for v in valid
        if v > 0
            return false
        end
    end
    return true
end

"""
    set_unconstrained!(model::Hildrethmodel)

Sets the model's solution element (Soln) to the unconstrained solution (ucSoln)
"""
function set_unconstrained!(model::HildrethModel)
    model.Soln = model.ucSoln
end

function is_empty(model::HildrethModel)

    E() = try
              return model.E
          catch UndefRefError
              return []
          end
    F() = try
              return model.F
          catch UndefRefError
              return []
          end

    return E() == [] && F() == [] &&
           model.constraints == OrderedDict() &&
           model.variable_count == 0 &&
           model.constraint_count == 0 
end

"""
    SC_HildrethConvergence(limit)

A stop condition that implements the following convergence check:
```math
(λ-λ_old)^T*(λ-λ_old) < limit
``` 
Where 'limit' is the value below which computation will halt. 
"""
struct SC_HildrethConvergence <: StoppingCondition
    value ::Float64
end

"""
    stopcondition(model::HildrethModel, convergence::SC_HildrethConvergence)

Convergence checking for original Hildreth implementation.
"""
function stopcondition(model::HildrethModel, 
                       convergence::SC_HildrethConvergence
                       )::Bool
    λ = model.workingvars["λ"]
    λ_old = model.workingvars["λ_old"]
    return (λ-λ_old)' * (λ-λ_old) < convergence.value, RAM_OPTIMAL
end

function get_iterations(model::HildrethModel)::Int
    return model.workingvars["iterations"]
end

"""
    resolver!(model::HildrethModel)

Calculates the primal result after convergence is met with the equation:

    x = -(E_f\\(F + M^T * λ))

for the problem: min 1/2(x'Ex) + F'x s.t. Mx <= γ

Where E_f is the factorised form of E in the problem statement, and λ is
a working variable of the row action method.
"""
function resolver!(model::HildrethModel)
    Ef = model.E_fact
    F = model.F
    Mt = get_constraintmatrix(model)
    λ = model.workingvars["λ"]
    model.Soln = -(Ef\(F + Mt * λ))
end

function answer(model::HildrethModel)
    if model.Soln == nothing 
        throw(ErrorException("Attempt to access answer value before any iterations have completed."))
    else
        return model.Soln
    end
end


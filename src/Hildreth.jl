import Base.==
export Hildreth, SC_HildrethConvergence

"""
    Hildreth(E, F, M, γ, H K, ucSoln, Soln, E_fact, workingvars, options)

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
mutable struct Hildreth{T} <: ModelFormulation
    #QP matrix
    E::Array{T}
    #QP vector
    F::Vector{T}
    H::Array{T}
    K::Vector{T}
    ucSoln::Vector{T}
    Soln::Vector{T}
    E_fact::Union{Bidiagonal,Factorization,Array,Diagonal}
    #TODO types
    λ
    λ_old

    function Hildreth{T}() where T
        model = new()
        return model
    end
end


"""
    GetModel(::Hildreth)::Hildreth
    
Returns a seleton model of the problem for solving with Hildreth's orignal
method.
"""
function GetModel(::Hildreth)::Hildreth
    return Hildreth()
end

"""
    iterate!(model::Hildreth)

Performs one iteration of the algorithm. Updates λ as it progresses.

Treats the entire summation as a calculation of H_i * λ, then subtracts the 
contribution of the currently considered λ. 
"""
function iterate!(model::Hildreth)
    λ = model.λ
    model.λ_old = copy(λ)
    H = model.H
    K = model.K
    
    for (i,l) in enumerate(λ)
        w = (H[i:i,:] * λ)[1] - H[i,i] * l
        w += K[i]
        w /= -H[i,i]
        λ[i] = max(0, w)
    end
end

"""
    setobjective!(model::Hildreth, E::Array{T, 2}, F::Vector{T}, num_vars::Int) where T

Sets the objective function for the problem. Hildreth's algorithm requires that 
E is positive definite. The number of problem variables is also needed.
"""
function setobjective!(model::Hildreth, E::Array{T, 2}, F::Vector{T}, num_vars::Int) where T
    model.E = E
    model.F = F
    model.status.variable_count = num_vars
end


function shrinkobjective(model::Hildreth, index::Int)
    model.E = model.E[setdiff(1:end, index), setdiff(1:end, index)]
    deleteat!(model.F, index)
end

"""
    buildmodel(model::Hildreth)

Builds the internal variables based on problem specification
"""
function buildmodel!(model::Hildreth)
    Mt = get_constraintmatrix(model)
    M = Mt'
    γ = get_constraintvector(model)

    model.E_fact    = cholesky(model.E)
    model.H         = M * (model.E_fact\Mt)
    model.K         = γ + (M * (model.E_fact\model.F))
    model.ucSoln    = -(model.E_fact\model.F)
    model.λ         = zeros(model.status.constraint_count)
    model.λ_old     = zeros(model.status.constraint_count)
end

"""
    ==(a::Hildreth, b::Hildreth)::Bool

Checks equality of two hildreth model structs. Does not check for object
equality, only value equality.
"""
#FIXME Update to  new model
function ==(a::Hildreth, b::Hildreth)::Bool
    return a.E == b.E &&
           a.F == b.F &&
           a.E_fact == b.E_fact &&
           a.H == b.H &&
           a.K == b.K &&
           a.ucSoln == b.ucSoln &&
           a.λ == b.λ &&
           a.λ_old == b.λ_old
end

"""
    valid_unconstrained(model::Hildreth)::Bool

Returns true if the initial minimum is within the constraints.
"""
function valid_unconstrained(model::Hildreth)::Bool
    valid = get_constraintmatrix(model)' * model.ucSoln - get_constraintvector(model)
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
function set_unconstrained!(model::Hildreth)
    model.Soln = model.ucSoln
end

function is_empty(model::Hildreth)
    #Best solution I could find quickly, probably better ways of doing it.
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

    @debug "Empty model test results:" E = E() F = F() status = empty_model_status(model)

    return E() == [] && F() == [] && empty_model_status(model)
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
    stopcondition(model::Hildreth, convergence::SC_HildrethConvergence)

Convergence checking for original Hildreth implementation.
"""
function stopcondition(model::Hildreth, 
                       convergence::SC_HildrethConvergence
                       )::Bool
    λ = model.λ
    λ_old = model.λ_old
    return (λ-λ_old)' * (λ-λ_old) < convergence.value, RAM_OPTIMAL
end

"""
    resolver!(model::Hildreth)

Calculates the primal result after convergence is met with the equation:

    x = -(E_f\\(F + M^T * λ))

for the problem: min 1/2(x'Ex) + F'x s.t. Mx <= γ

Where E_f is the factorised form of E in the problem statement, and λ is
a working variable of the row action method.
"""
function resolver!(model::Hildreth)
    Ef = model.E_fact
    F = model.F
    Mt = get_constraintmatrix(model)
    λ = model.λ
    model.Soln = -(Ef\(F + Mt * λ))
end

function variable_values(model::Hildreth)
    if model.Soln == nothing 
        throw(ErrorException("Attempt to access answer value before any iterations have completed."))
    else
        return model.Soln
    end
end

function objective_value(model::Hildreth)
    if model.Soln == nothing 
        throw(ErrorException("Attempt to access answer value before any iterations have completed."))
    else
        return 0.5 * (model.Soln' * model.E * model.Soln) + (model.Soln' * model.F)
    end
end


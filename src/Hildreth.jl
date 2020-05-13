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
mutable struct Hildreth{T} <: ModelFormulation
    #QP matrix
    E::SparseMatrixIndex
    #QP vector
    F::SparseVectorIndex

    H::SparseMatrixIndex
    K::SparseVectorIndex
    ucSoln::SparseVectorIndex
    Soln::SparseVectorIndex
    E_fact::SparseMatrixIndex
    λ::SparseVectorIndex
    λ_old::SparseVectorIndex

    Hildreth{T}() where T = new()

end

ObjectiveType(::Hildreth) = Quadratic()

"""
    iterate!(model::Hildreth)

Performs one iteration of the algorithm. Updates λ as it progresses.

Treats the entire summation as a calculation of H_i * λ, then subtracts the 
contribution of the currently considered λ. 
"""
iterate_args(::Hildreth) = [:λ, :λ_old, :H, :K]

function Iterate(method::Hildreth, model::RAMProblem)
    λ = GetSparse(model.λ)

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
    Setup

"""
function Setup(method::Hildreth, E, F)
    method.E = E
    method.F = F
end


#TODO ref update
function shrinkobjective(model::Hildreth, index::Int)
    model.E = model.E[setdiff(1:end, index), setdiff(1:end, index)]
    deleteat!(model.F, index)
end

"""
    buildmodel(model::Hildreth)

Builds the internal variables based on problem specification
"""
function Build(method::Hildreth, model)
    Mt = get_constraintmatrix(model)
    M = Mt'
    γ = get_constraintvector(model)
    E = GetSparse(model, method.E)
    F = GetSparse(model, method.F)

    method.E_fact    = RegisterSparse(model, cholesky(E))
    E_fact           = GetSparse(model, method.E_fact)

    method.H         = RegisterSparse(model, M * (E_fact\Mt))
    method.K         = RegisterSparse(model, γ + (M * (E_fact\F)))
    method.ucSoln    = RegisterSparse(model, -E_fact\F)

    method.λ         = RegisterSparse(model, zeros(model.constraint_count))
    method.λ_old     = RegisterSparse(model, zeros(model.constraint_count))
end

"""
    valid_unconstrained(model::Hildreth)::Bool

Returns true if the initial minimum is within the constraints.
"""
function valid_unconstrained(model::Hildreth)::Bool
    valid = get_constraintmatrix(model)' * GetSparse(model.ucSoln) - get_constraintvector(model)
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

    return GetSparse(E()) == [] && GetSparse(F()) == [] && empty_model_status(model)
end

"""
    SC_HildrethConvergence(limit)

A stop condition that implements the following convergence check:
```math
(λ-λ_old)^T*(λ-λ_old) < limit
``` 
Where 'limit' is the value below which computation will halt. 
"""
#TODO naming convention
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
    Ef = GetSparse(model.E_fact)
    F  = GetSparse(model.F)
    Mt = GetSparse(get_constraintmatrix(model))
    λ  = GetSparse(model.λ)
    model.Soln = RegisterSparse(-(Ef\(F + Mt * λ)))
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


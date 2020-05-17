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
    λ::Vector{T}
    λ_old::Vector{T}

    Hildreth{T}() where T = new()

end

ObjectiveType(::Hildreth) = Quadratic()
SupportsVariableDeletion(::Hildreth) = true 

function DeleteVariable(method::Hildreth, index)
    method.λ[index] = 0.0
    method.λ_old[index] = 0.0
end

"""
    Iterate(model::Hildreth)

Performs one iteration of the algorithm. Updates λ as it progresses.

Treats the entire summation as a calculation of H_i * λ, then subtracts the 
contribution of the currently considered λ. 
"""
iterate_args(::Hildreth) = [:H, :K]
function Iterate(model::RAMProblem, method::Hildreth, H, K)
    method.λ_old = method.λ
    for i in 1:model.constraint_count
        w = (H[i:i,:] * method.λ)[1] - H[i,i] * method.λ[i]
        w += K[i]
        w /= -H[i,i]
        method.λ[i] = max(0, w)
    end
end

"""
    Setup

"""
function Setup(method::Hildreth, E, F)
    method.E = E
    method.F = F
end

"""
    Build(model::Hildreth)

Builds the internal variables based on problem specification
"""
function Build(method::Hildreth, model)
    Mt = get_constraintmatrix(model)
    M = Mt'
    γ = get_constraintvector(model)
    E = GetSparse(model, method.E)
    F = GetSparse(model, method.F)
    
    #TODO: Likely an iterative solution that avoids casting F to a dense matrix
    F_v = Vector(F)
    #TODO re-implement E_fact
    #method.E_fact    = RegisterSparse(model, cholesky(E))
    #E_fact           = GetSparse(model, method.E_fact)

    method.H         = RegisterSparse(model, M * (E\Mt))
    method.K         = RegisterSparse(model, γ + (M * (E\F_v)))
    method.ucSoln    = RegisterSparse(model, -E\F_v)

    method.λ         = zeros(model.constraint_count)
    method.λ_old     = zeros(model.constraint_count)
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

function is_empty(method::RAMProblem, model::Hildreth)
    #Best solution I could find quickly, probably better ways of doing it.
    entry(sym) = try
                    return GetSparse(model, getproperty(method, sym))
                 catch UndefRefError
                    return []
                 end

    return isempty(entry(:E)) && isempty(entry(:F))
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
resolve_args(::Hildreth) = [:E, :F]
function Resolve(model::RAMProblem, method::Hildreth, E, F)
    Mt = get_constraintmatrix(model)
    method.Soln = RegisterSparse(model, -(E\Vector(F + Mt * method.λ)))
end

GetVariables(model::RAMProblem, method::Hildreth) = GetSparse(model, method.Soln)

function objective_value(model::Hildreth)
    if model.Soln == nothing 
        throw(ErrorException("Attempt to access answer value before any iterations have completed."))
    else
        return 0.5 * (model.Soln' * model.E * model.Soln) + (model.Soln' * model.F)
    end
end


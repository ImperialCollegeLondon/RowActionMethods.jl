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
mutable struct HildrethModel <: ModelFormulation
    E::Array{Float64}
    F::Vector{Float64}
    M::Array{Float64}
    γ::Vector{Float64}
    H::Array{Float64}
    K::Vector{Float64}
    ucSoln::Vector{Float64}
    Soln::Union{Vector{Float64},Nothing}
    E_fact::Union{Factorization,Array,Diagonal,Nothing}
    workingvars::WorkingVars
end


"""
    GetModel(::Hildreth)::HildrethModel
    
Returns a seleton model of the problem for solving with Hildreth's orignal
method.
"""
function GetModel(::Hildreth)::HildrethModel
    return HildrethModel([],[],[],[],[],[],[],
                         Nothing(),Nothing(), 
                         Dict("iterations"=>0))
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
    buildmodel(E, F, M, γ, ::Hildreth)

Returns the problem model for the original Hildreth method.
"""
function buildmodel!(model::HildrethModel, E, F, M, γ)
    model.E = E
    model.F = F
    model.M = M
    model.γ = γ
    model.E_fact = factorize(E)
    model.H = M * (model.E_fact\transpose(M))
    model.K = γ + (M * (model.E_fact\F))
    model.ucSoln = -(model.E_fact\F)
    model.workingvars["λ"] = zeros(size(model.γ))
    model.workingvars["λ_old"] = zeros(size(model.γ))
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
    return model.E == [] &&
           model.F == [] &&
           model.M == [] &&
           model.γ == []
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
    stopcondition(vars::WorkingVars, convergence::SC_HildrethConvergence)

Convergence checking for original Hildreth implementation.
"""
function stopcondition(model::HildrethModel, 
                       convergence::SC_HildrethConvergence
                       )::Bool
    λ = model.workingvars["λ"]
    λ_old = model.workingvars["λ_old"]
    return (λ-λ_old)' * (λ-λ_old) < convergence.value
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
    M = model.M
    λ = model.workingvars["λ"]
    model.Soln = -(Ef\(F + M' * λ))
end

function answer(model::HildrethModel)
    if model.Soln == nothing 
        throw(ErrorException("Attempt to access answer value before any iterations have completed."))
    else
        return model.Soln
    end
end

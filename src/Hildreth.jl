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
    A::SparseMatrixCSC{T}
    b::SparseVector{T}

    Δ::SparseMatrixCSC{T}
    n::Int

    z::SparseVector{T}
    x::SparseVector{T}

    user_initial_point::Union{Nothing,Vector{T}}

    function Hildreth{T}(;kwargs...) where T
        m = new()
        m.user_initial_point = 
            haskey(kwargs, :initial) ? kwargs[:initial] : nothing
        return m
    end
end

ObjectiveType(::Hildreth) = Quadratic()
SupportsVariableDeletion(::Hildreth) = true 

function DeleteVariable(method::Hildreth, index)
    method.λ[index] = 0.0
    method.λ_old[index] = 0.0
end

"""
    Build(model::Hildreth)

Builds the internal variables based on problem specification
"""
function Build(model::RAMProblem, method::Hildreth)
    G = get_constraintmatrix(model)'
    h = get_constraintvector(model)

    B, d = GetObjectiveFactorised(model)
    
    method.A        = G/B.U
    method.b        = h + (G/B)d

    method.Δ        = (method.A * method.A')'

    method.n        = length(h)

    #TODO adjust this range depending on input value
    if method.user_initial_point == nothing
        method.z = rand(0.1:0.1:10.0, method.n)
    else
        method.z = method.user_initial_point
    end
    method.x        = -method.A'method.z
end

"""
    Iterate(model::Hildreth)

Performs one iteration of the algorithm. Updates λ as it progresses.

Treats the entire summation as a calculation of H_i * λ, then subtracts the 
contribution of the currently considered λ. 
"""
iterate_args(::Hildreth) = [:z, :Δ, :b, :n]
function Iterate(model::RAMProblem, method::Hildreth, z, Δ, b, n)
    for i in 1:n
        w = Δ[:,i]'*z
        w += b[i]
        w /= Δ[i,i]
        w = z[i] - w
        method.z[i] = max(0,w)
    end
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

Where E_f is the factorised form of E in the problem statement, and λ is
a working variable of the row action method.
"""
resolve_args(::Hildreth) = [:A, :z]
function Resolve(model::RAMProblem, method::Hildreth, A, z)
    method.x = -A'*z
end

function GetVariables(model::RAMProblem, method::Hildreth) 
    B, d = GetObjectiveFactorised(model)
    return B.U\method.x - B\Vector(d)
end


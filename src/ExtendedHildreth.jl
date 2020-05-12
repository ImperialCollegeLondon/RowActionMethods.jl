using LinearAlgebra
export ExtendedHildreth

#TODO Implement almost-cyclic functionality
#TODO Vectorise implementation of algorithm - possibly unneeded?
#TODO Find neater method of initialisaiton than unions of desired types and nothing
#TODO Add methods to check convergence.

"""
    ExtendedHildreth

A model for solving the problem with the extended Hildreth's method,
presented by Lent and Censor.

Problem in the form:
min ½⟨By,y⟩+⟨y,d⟩ s.t. Gy≦h

Reformulated to:
min ½∥x∥² s.t. Ax≦b

where:
A=GD^-1
b=h+GB^-1

from:
B=D'D (Cholesky decomposition)
y=inv(D)x-inv(B)d
"""
mutable struct ExtendedHildreth <: ModelFormulation
    #Original problem
    y::Union{Vector{Float64},Nothing}#Solution
    B::Array{Float64}
    d::Vector{Float64}
    G::Array{Float64}
    h::Vector{Float64}
    n::Int32
    m::Int32
    #Reformulation
    B_chol::Union{Cholesky{Float64},Nothing}#Defined as union for initialisation
    D::Union{UpperTriangular{Float64},Nothing}#Defined as union for initialisation
    A::Array{Float64}
    b::Vector{Float64}
    #Working vars 
    iterations::Int
    x::Vector{Float64}
    z::Vector{Float64}
    #Parameters
    relaxation::Function
    z_initial::Union{Vector{Number}, Number}
end
"""
    GetModel(::ExtendedHildreth, z_initial::Union{Vector{Number}, Number}, 
    relaxation_func::Function)::ExtendedHildreth

Returns a skeleton model of the problem for solving with the extended Hildreth's modethod. 

z_initial is the initial value of the z dual variable. Any supplied vector should be of the correct dimensions (z∈ℝᵐ, where m is the number of constraints). If a scalar k is supplied then z will be a column vector of ones scaled by k.

This method supports the inclusion of a relaxation parameter, this is a scaling that is applied on each iteration (see the original paper for details). This passed function should calculate the relaxation parameter value from the current iteration index. All resultant values should be close to 1 for convergence to hold. By default relaxation will be 1.
"""
function GetModel(::ExtendedHildreth, 
                   z_initial::Union{Vector{Number}, Number},
                   relaxation_func::Function
                  )::ExtendedHildreth
    return ExtendedHildreth(Nothing(),[],[],[],[],0,0,
                                 Nothing(),Nothing(),[],[],
                                 0,[],[],
                                 relaxation_func, z_initial)
end

#Better ideas on how to structure this are welcome 
"""
    GetModel(::ExtendedHildreth; options...)

Optional arguments:
z_initial - initial value for z to take, must be in 
the non-negative orthant of Rⁿ (defaults to vector of ones). Can be 
a vector (of the correct dimension for the problem) or a scalar multiplier of the vector of ones.
relaxation - a value for the relaxation series to take, either a constant or a function to generate relaxation value from the iteration index. All values should be positive for convergence proof to hold.
"""
function GetModel(::ExtendedHildreth; options...)
    options = Dict(options)
    if haskey(options, "z_initial") 
        z=options["z_initial"]
    else
        z=1
    end

    if haskey(options, "relaxation")
        r = options["relaxation"]
        if typeof(r) == Number
            relaxation = (k->r)
        elseif typeof(r) == Function
            relaxation = r
        else
            error("Invalid relaxation type, relaxation should be passed as a function or Number")
        end
    else
        relaxation = (k->1)
    end

    return GetModel(ExtendedHildreth(), z, relaxation)
end

function buildmodel!(model::ExtendedHildreth, B, d, G, h)
    model.B=B
    model.d=d
    model.G=G
    model.h=h

    model.n=size(d)[1]
    model.m=size(h)[1]
    
    model.B_chol=cholesky(B)
    model.D=model.B_chol.U
    model.A=G/model.D #A=GD^-1
    model.b=h+(G * (model.B_chol\d)) #b = h + GB^-1d

    #If z_intial is a number, give vector of number
    if typeof(model.z_initial) <: Number
        model.z=model.z_initial * zeros(model.m)
    #Error if z_initial is a vector of the wrong dimension
    elseif size(model.z_initial) != model.m
        error("Problem definition does not fit supplied initial value for z")
    else 
        model.z = model.z_inital
    end

    model.x=-model.A'*model.z 
end

#slightly unsure about the algorithm 
function iterate!(model::ExtendedHildreth)
    #c = min(zᵢ, uᵢ)
    #uᵢ = r * (bᵢ-⟨aᵢ,x⟩)/∥aᵢ∥²
    #xₖ₊₁ = xₖ+c_aᵢ
    #zₖ₊₁ = zₖ+c_eᵢ
    
    r = model.relaxation(model.iterations)
    b = model.b
    A = model.A
    x = model.x
    z = model.z
    n = model.n
    m = model.m

    x_temp = zeros(n)
    z_temp = zeros(m)

    for i = 1:m
        u = r * (b[i] - A[i,:]⋅x) / norm(A[i,:])^2
        c = min(z[i], u)
        x_temp += (c .* A[i,:])
        z_temp[i] = c
    end

    model.x += x_temp
    model.z -= z_temp 
    model.iterations += 1
    
end

"""
    validmodel(model::ExtendedHildreth)::Bool

Method doesn't implement a prior check for viability, therefore function
always returns false to ensure that solver runs. 
"""
function valid_unconstrained(model::ExtendedHildreth)::Bool
    return false
end

"""
    get_unconstrained(model::ExtendedHildreth)

As a prior check for viability is never resolved as true, this method 
does nothing.
"""
function get_unconstrained(model::ExtendedHildreth)
end

function get_iterations(model::ExtendedHildreth)::Int
    return model.iterations
end

function resolver!(model::ExtendedHildreth)
    model.y = model.D\model.x - model.B_chol\model.d
end

function answer(model::ExtendedHildreth)
    if model.y == nothing
        throw(ErrorException("Attempt to access answer value before any iterations have completed."))
    else
        return model.y
    end
end


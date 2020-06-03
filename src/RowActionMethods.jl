module RowActionMethods

using DataStructures
using SparseArrays
using Base.Threads: @threads, nthreads

#=
Author: Edward Stables 
Date: 26-11-2019

For solving QP problems of the form:
Minimize: J = 1/2(x'Ex) + F'x
Subject to the constraints: Mx <= γ
=#

using LinearAlgebra
using Logging

const RAM = RowActionMethods
export RAM

include("./Types.jl")
include("./StopConditions.jl")
include("./MOI_wrapper.jl")
include("./Constraints.jl")
include("./Core.jl")

#Method includes 
include("./Hildreth.jl")
include("./ExtendedHildreth.jl")

#Maps strings to method types
method_mapping= Dict("Hildreth" => Hildreth,
                     "ExtendedHildreth" => ExtendedHildreth,
                    )
end # module

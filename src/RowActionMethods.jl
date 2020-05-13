module RowActionMethods

using DataStructures
using SparseArrays

#TODO Standardise notation between methods
#=
Author: Edward Stables 
Date: 26-11-2019

For solving QP problems of the form:
Minimize: J = 1/2(x'Ex) + F'x
Subject to the constraints: Mx <= γ
=#

using LinearAlgebra
using Logging

export iterate_model!, GetModel, buildmodel!, answer, get_SC

const RAM = RowActionMethods
export RAM

include("./Types.jl")
include("./StopConditions.jl")
include("./MOI_wrapper.jl")
include("./Benchmarks.jl")
include("./Constraints_Variables.jl")
include("./Core.jl")

#Method includes 
include("./Hildreth.jl")
include("./ExtendedHildreth.jl")

#Maps strings to method types
method_mapping= Dict("Hildreth" => Hildreth,
                     "ExtendedHildreth" => ExtendedHildreth,
                    )
function run_benchmarks()
    run_all_benchmarks()
end

end # module

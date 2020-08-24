module RowActionMethods

using SuiteSparse
using DataStructures
using SparseArrays
using Base.Threads: @threads, nthreads
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

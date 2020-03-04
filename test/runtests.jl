using RowActionMethods
using LinearAlgebra
using Test

const RAM = RowActionMethods

include("example_qp.jl")

include("test_API.jl")
#include("test_MOI_wrapper.jl")
include("test_Hildreth.jl")

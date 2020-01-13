using MathOptInterface
const MOI = MathOptInterface
const MOIT = MOI.Test
const MOIU = MOI.Utilities
const MOIB = MOI.Bridges

#Test for single method as the interface should apply to all, if other tests pass
const optimizer = RAM.Optimizer("Hildreth")
MOI.set(optimizer, MOI.Silent(), true)

@testset "SolverName" begin
    @test MOI.get(optimizer, MOI.SolverName()) == "RowActionMethods-Hildreth"
end

@testset "moi supports_default_copy_to" begin
    @test MOIU.supports_default_copy_to(optimizer, false)
    # Use `@test !...` if names are not supported
    @test MOIU.supports_default_copy_to(optimizer, true)
end

const bridged = MOIB.full_bridge_optimizer(optimizer, Float64)
const config = MOIT.TestConfig(atol=1e-6, rtol=1e-6)

#=
@testset "MOI Unit" begin
    MOIT.unittest(bridged, config)
end

@testset "MOI Modification" begin
    MOIT.modificationtest(bridged, config)
end
=#

@testset "MOI Continuous Linear" begin
    MOIT.contlineartest(bridged, config)
end


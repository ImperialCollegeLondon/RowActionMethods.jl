@testset "Core" begin

    @testset "GetModel" begin 
        #Only test for basic functionality here.
        #Put more detailed tests in individual solver test files.

        @test_throws MethodError RAM.GetModel()
        @test_throws MethodError RAM.GetModel(;some_option = true)
 
        #please don't name a solver this
        @test_throws ArgumentError RAM.GetModel("Nonexistent Algorithm") 

        problem = RAM.GetModel("Hildreth")
        @test typeof(problem) == RAM.RAMProblem{Float64,Int64}
        @test typeof(problem.method) == RAM.Hildreth{Float64}

    end
end

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

    @testset "_objective_type" begin
        problem = RAM.GetModel("Hildreth")
        @test RAM._objective_type(problem) == RAM._objective_type(problem.method) == RAM.Quadratic()
    end

    @testset "SetThreads" begin
        problem = RAM.GetModel("Hildreth")
        @test problem.threads == false

        RAM.SetThreads(problem)
        @test problem.threads == true

        RAM.SetThreads(problem; threads=true)
        @test problem.threads == true

        RAM.SetThreads(problem; threads=false)
        @test problem.threads == false
    end

    @testset "GetObjective" begin

    end

    @testset "GetObjectiveFactorised" begin end

    @testset "Iterate" begin end

    @testset "Resolve" begin end
    @testset "GetVariables" begin end
    @testset "is_empty" begin end
    @testset "get_model_status" begin end
    @testset "objective_value" begin end
    @testset "SupportsDeleteConstraint" begin end
    @testset "IterateRow" begin end
    @testset "GetTempVar" begin end
    @testset "VarUpdate" begin end
    @testset "objective_value" begin end

    @testset "Setup" begin end
    @testset "_init_run" begin end
    @testset "optimize!" begin end


end

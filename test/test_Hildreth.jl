@testset "Hildreth" begin
    #Standard Tests
    RowActionMethodStandardTests(Hildreth())

    empty_model = RowActionMethods.HildrethModel([],[],[],[],[],[],[],
                                                 Nothing(),Nothing(),
                                                 Dict("iterations"=>0))
    built_model = RowActionMethods.HildrethModel(
                  [1.0 0.0; 0.0 1.0], [-2.0, -2.0], [1.0 0.0; -1.0 0.0; 0.0 1.0; 0.0 -1.0], 
                  [1.0, 0.0, 1.0, 0.0], 
                  [1.0 -1.0 0.0 0.0; -1.0 1.0 0.0 0.0; 0.0 0.0 1.0 -1.0; 0.0 0.0 -1.0 1.0], 
                  [-1.0, 2.0, - 1.0, 2.0], [2.0, 2.0], nothing, [1 0; 0 1], 
                  Dict{String,Any}("iterations" => 0,"λ" => [0.0, 0.0, 0.0, 0.0], 
                                   "λ_old" => [0.0, 0.0, 0.0, 0.0])
                 )
    
    problem = example_2()
    question = problem["prob"]
    
    @testset "Model" begin
        model = Optimizer(Hildreth()) 
        @test model == empty_model
        buildmodel!(model, question...)
        @test model == built_model
        @test empty_model != built_model
    end

    @testset "Iterations" begin
        #TODO: Implement this with some example problems that need more iterations
    end
    
    @testset "Functions" begin
        built_model.ucSoln = [1; 1]    
        @test RowActionMethods.valid_unconstrained(built_model) == true
        built_model.ucSoln = [2; 2]    
        @test RowActionMethods.valid_unconstrained(built_model) == false

        RowActionMethods.set_unconstrained!(built_model)
        @test built_model.Soln == built_model.ucSoln

        @test answer(built_model) == built_model.Soln
    end

end

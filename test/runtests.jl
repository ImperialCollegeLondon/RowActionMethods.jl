using RowActionMethods
using LinearAlgebra
using Test


"""
    RowActionMethodStandardTests(method::T,
                                 run_answer_check=true
                                ) where T<:RowActionMethod

Contains a series of standard testsets that each row action method should pass 
if following the intended design. This should be run in addition to any custom
tests needed by a method, and is intended to ensure compatibility and check 
for common errors.

This mostly focuses on the behaviour in regards to the iteration stopping 
conditions. It also checks the results of the optimisation of a simple problem
after 20 iterations. 20 iterations should be sufficient for solving, but if
this poses problems for a method then disable by passing false to
run_answer_check. If this is done then please ensure your tests implement a 
similar check on a simple problem with the needed setup.
"""
function RowActionMethodStandardTests(method::T; 
                                      run_answer_check=true
                                     ) where T<:RowActionMethods.RowActionMethod
    testset_name = string(typeof(method))
    #TODO move these example problems to tests 
    question = RowActionMethods.example_1()
    p = question["prob"]
    a = question["ans"]
    
    @testset "$testset_name" begin
        model = Optimizer(method)
        @testset "Optimizer" begin
            @test typeof(model) <: RowActionMethods.ModelFormulation
        end

        buildmodel!(model, p...)
        @testset "BuildModel" begin
            @test RowActionMethods.get_iterations(model) == 0
            #Suggestions on a more appropriate error are welcome
            @test_throws ErrorException answer(model)
        end

        @testset "Single Iteration" begin
             RowActionMethods.iterate!(model)
             @test RowActionMethods.get_iterations(model) == 1
        end


        @testset "Iterations Stop Condition" begin
            iterations_condition = SC_Iterations(11)  
            iterate_model!(model, iterations_condition)
            @test RowActionMethods.get_iterations(model) == 11

            iterations_condition = SC_Iterations(5)
            iterate_model!(model, iterations_condition)
            @test RowActionMethods.get_iterations(model) == 11

            iterations_condition = SC_Iterations(20)
            iterate_model!(model, iterations_condition)
            @test RowActionMethods.get_iterations(model) == 20
        end

       if run_answer_check
            @testset "Answer Check" begin
                @test answer(model) ≈ a
            end
       end
    end
end

RowActionMethodStandardTests(Hildreth())
RowActionMethodStandardTests(ExtendedHildreth())

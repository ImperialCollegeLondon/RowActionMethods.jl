using RowActionMethods
using LinearAlgebra
using Test

const RAM = RowActionMethods

include("example_qp.jl")

@testset "RowActionCore" begin
    @testset "get_SC" begin
        a = SC_Iterations(5)
        b = SC_Iterations(10)
        c = SC_Iterations(15)

        result_1 = MultipleStopCondition(RAM.StoppingCondition[
                                         SC_Iterations(5)
                                         ])
        result_2 = MultipleStopCondition(RAM.StoppingCondition[
                                         SC_Iterations(5),
                                         SC_Iterations(10)
                                         ])
        result_3 = MultipleStopCondition(RAM.StoppingCondition[
                                         SC_Iterations(5),
                                         SC_Iterations(10),
                                         SC_Iterations(15)
                                         ])

        @test RAM.compare_MultipleStopCondition(get_SC(a), result_1)
        @test RAM.compare_MultipleStopCondition(get_SC(a, b), result_2)
        @test RAM.compare_MultipleStopCondition(get_SC(a, b, c), result_3)
    end
end

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

TODO: Add tests to check for behaviour with different kinds of input matrices,
ideally should test every type that the `factorize()` function returns.  
"""
function RowActionMethodStandardTests(method::T; 
                                      run_answer_check=true
                                     ) where T<:RAM.RowActionMethod
    testset_name = string(typeof(method), "_Standard")
    #TODO move these example problems to tests 
    question = example_1()
    p = question["prob"]
    a = question["ans"]
    
    @testset "$testset_name" begin
        model = GetModel(method)
        @testset "GetModel" begin
            #Check model is correct type
            @test typeof(model) <: RAM.ModelFormulation
        end

        buildmodel!(model, p...)
        @testset "BuildModel" begin
            #Check model is initialised with 0 iterations
            @test RAM.get_iterations(model) == 0
            #Check that an answer isn't present if no iterations have run
            #Suggestions on a more appropriate error are welcome
            @test_throws ErrorException answer(model)
        end

        @testset "Single Iteration" begin
            RAM.iterate!(model)
            #Check an iteration updates common working variables
            @test RAM.get_iterations(model) == 1
        end


        @testset "Iterations Stop Condition" begin
            #Check that iteration stops on a given limit
            iterations_condition = SC_Iterations(11)  
            iterate_model!(model, iterations_condition)
            @test RAM.get_iterations(model) == 11

            #Check that no iterations are performed if the limit is reached
            iterations_condition = SC_Iterations(5)
            iterate_model!(model, iterations_condition)
            @test RAM.get_iterations(model) == 11

            #Check that more iterations are performed on a limit update
            iterations_condition = SC_Iterations(20)
            iterate_model!(model, iterations_condition)
            @test RAM.get_iterations(model) == 20
        end

       if run_answer_check
            @testset "Answer Check" begin
                @test answer(model) ≈ a
            end
       end
    end
end

include("test_MOI_wrapper.jl")
#include("test_Hildreth.jl")
#include("test_ExtendedHildreth.jl")

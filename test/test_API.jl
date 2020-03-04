#TODO implement function when ModelFormulation supports parametric types
#get_model_parameter(::RAM.ModelFormulation{T}) where {T} = T
get_status_parameter(::RAM.RAM_Components{T}) where {T} = T

"""
    test_API(method::RowActionMethod)

Checks each required function within the API for a given solver. 

"""
#TODO: Add option to test problem modification functionality 
function test_API(method::RAM.RowActionMethod)
    testset_name = string(typeof(method), "_api_test")

    @testset "$testset_name" begin

        model = RAM.GetModel(method)
        @testset "ModelFormulation" begin
            @test typeof(model) <: RAM.ModelFormulation
            
            @test hasfield(typeof(model), :status)
            param = get_status_parameter(model.status)
            @test typeof(model.status) == RAM.RAM_Components{param}
            @test RAM.is_empty(model)
        end

        E = [20.0 6.0; 6.0 5.0]
        F = [3.0; 2.0]
        Q = [1.0;2.0]
        b = 4.0
        N = 2
        ans = [ -0.046875; -0.34375]
        cost = -0.4140625

        @testset "set_objective!" begin
            model = RAM.GetModel(method)
            RAM.setobjective!(model, E, F, N)
            @test model.status.variable_count == N
        end
        
        con = RAM.addconstraint!(model, Q, b) 
        buildmodel!(model)
        iterate_model!(model)        

        @testset "return values" begin
            @test RAM.variable_values(model) ≈ ans
            @test RAM.objective_value(model) ≈ cost
        end
    end
end

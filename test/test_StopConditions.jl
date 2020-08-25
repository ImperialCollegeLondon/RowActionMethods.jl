@testset "Stop Conditions" begin

    @testset "Iteration Limit" begin
        @test_throws DomainError IterationCondition(0)
        @test_throws DomainError IterationCondition(-1)

        ic    = IterationCondition(3)
        model = RAM.RAMProblem{Float64, Int64}("Hildreth")
        model.iterations = 0

        # Verify the condition works
        @test !RAM.test_stopcondition(model, ic)

        model.iterations = 1
        @test !RAM.test_stopcondition(model, ic)

        model.iterations = 2
        @test !RAM.test_stopcondition(model, ic)

        model.iterations = 3
        @test RAM.test_stopcondition(model, ic)

        # Verify it has the right status code
        @test RAM.status_code(ic) == RAM.ITERATION_LIMIT_REACHED
    end

    @testset "Time Limit" begin
        @test_throws DomainError TimeCondition(0.0)
        @test_throws DomainError TimeCondition(-1.0)

        tc    = TimeCondition(5.0)
        model = RAM.RAMProblem{Float64, Int64}("Hildreth")

        # Verify it initialized properly
        @test tc.start ≈ 0.0
        RAM.init_stopcondition(model, tc)
        @test tc.start > 0.0

        # Verify the condition works
        @test !RAM.test_stopcondition(model, tc)

        sleep(2.0)
        @test !RAM.test_stopcondition(model, tc)

        sleep(4.0)
        @test RAM.test_stopcondition(model, tc)

        # Verify it has the right status code
        @test RAM.status_code(tc) == RAM.TIME_LIMIT_REACHED
    end

    @testset "Multiple conditions" begin
        let model = RAM.RAMProblem{Float64, Int64}("Hildreth"),
            stopVec = [ IterationCondition(3);
                        TimeCondition(5.0) ]

            # Ensure the TimeCondition started properly
            RAM._init_stopconditions(model, stopVec)
            @test stopVec[2].start > 0

            # Test that the first stop condition fires when it is supposed to
            model.iterations = 0
            @test RAM._check_stopconditions(model, stopVec)
            sleep(6.0)
            @test !RAM._check_stopconditions(model, stopVec)

            @test model.status == RAM.TIME_LIMIT_REACHED
        end

        # Rebuild the parts to try the other condition this time
        let model = RAM.RAMProblem{Float64, Int64}("Hildreth"),
            stopVec = [ IterationCondition(3);
                        TimeCondition(500.0) ]

            RAM._init_stopconditions(model, stopVec)

            model.iterations = 0
            @test RAM._check_stopconditions(model, stopVec)
            model.iterations = 1
            @test RAM._check_stopconditions(model, stopVec)
            model.iterations = 3
            @test !RAM._check_stopconditions(model, stopVec)
            @test model.status == RAM.ITERATION_LIMIT_REACHED
        end
    end

end

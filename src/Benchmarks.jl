using LinearAlgebra
using BenchmarkTools
using MathOptInterface

#include("ProblemList.jl")

#Currently only using the problems available in OptimizationProblems repo for 
#simplicity CUTEst naming convention is available here: 
#http://www.cuter.rl.ac.uk//Problems/classification.shtml
#This solver requires quadratic problems with linear constraints, therefore 
#classification of QLR{0,1,2}-.N-dd is suitable. This gives quite small problems
#from this problem set, larger are needed in future.

function run_all_benchmarks()
    single_benchmark(QLR2_AN_4_6)
end

function single_benchmark(problem)
    prob, var = problem()
    #Build problem build benchmark
    @benchmark set_optimizer($prob, with_optimizer(RAM.Optimizer, "Hildreth"))
    @benchmark optimize($prob)
end

"""
    function generate_random_mat(n::Int, ::Type{T}; entry_gen=rand, non_zero_gen=0.1)::Array{T,2} where T
    
Returns nxn matrix with non-zero entry selection defined by `entry_gen`, 
and probabilty of an entry being non-zero given by `non_zero_gen`. The first
entry in entry_gen should denote the return type.

"""
function generate_random_posdef_matrix(n::Int, ::Type{T}; entry_gen=rand, non_zero_gen::Float64=0.1)::Array{T,2} where T
    #FIXME: This function is terrible and needs to change
    M = zeros(T,n,n)

    while true
        for i=1:n,j=1:n
            if rand() < non_zero_gen
                M[i,j] = entry_gen(T)
            end
        end
        M = M'*M
        if isposdef(M)
            break
        end
        M = zeros(T, n, n)
    end

    return M
end

function generate_random_vector(n::Int, ::Type{T}; entry_gen=rand, non_zero_gen::Float64=0.1)::Vector{T} where T
    V = zeros(T,n)

    for i=1:n
        if rand() < non_zero_gen
            V[i] = entry_gen(T)
        end
    end

    return V
end


function moi_generate_test(model::Optimizer, n::Int, m::Int
                          )::Tuple{Vector{MOI.VariableIndex},Vector{MOI.ConstraintIndex}}
    #TODO make more generic to other solvers
    #TODO give specific value ranges (currently mostl limited to [0,1])
    
    x = MOI.add_variables(model, n)
    obj_quad_term = MOI.ScalarQuadraticTerm{Float64}[]
    obj_aff_term = MOI.ScalarAffineTerm{Float64}[]

    E = generate_random_posdef_matrix(n, Float64)
    F = generate_random_vector(n, Float64)

    for i=1:n, j=1:n
        push!(obj_quad_term, MOI.ScalarQuadraticTerm(E[i,j], x[i], x[j]))
    end

    for i=1:n
        push!(obj_aff_term, MOI.ScalarAffineTerm(F[i], x[i]))
    end

    obj = MOI.ScalarQuadraticFunction(obj_aff_term, obj_quad_term, 0.0)

    MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}}(), obj) 

    #Constraint generation
    con = MOI.ConstraintIndex[]
    for i = 1:m
        F = 10 .* generate_random_vector(n, Float64, non_zero_gen=0.5)
        func = MOI.ScalarAffineFunction(
               map((c, i) -> MOI.ScalarAffineTerm(c, i), F, x), 0.0)
        lim = MOI.LessThan(rand())
        push!(con, MOI.add_constraint(model, func, lim))
    end

    return x, con
end


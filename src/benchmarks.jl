#using BenchmarkTools
using JuMP

include("ProblemList.jl")

#Currently only using the problems available in OptimizationProblems repo for 
#simplicity CUTEst naming convention is available here: 
#http://www.cuter.rl.ac.uk//Problems/classification.shtml
#This solver requires quadratic problems with linear constraints, therefore 
#classification of QLR{0,1,2}-.N-dd is suitable. This gives quite small problems
#from this problem set, larger are needed in future.

const suitable_probs = [QLR2_AN_4_6]

function run_all_benchmarks()
    single_benchmark(QLR2_AN_4_6)
end

function single_benchmark(problem)
    prob, x = problem()

    set_optimizer(prob, with_optimizer(Optimizer, "Hildreth"))
    optimize!(prob)
    for v in x
        println(value(v))
    end
end

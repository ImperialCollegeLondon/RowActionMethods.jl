"""
Provides the input values to replicate example 2.10 (from the book refered to 
in the introductory comments)
"""
function example_1()::Dict
    E = [2 -1; -1 1]
    F = [-1; 0]
    M = [-1 0; 0 -1; 3 2]
    γ = [0;0;4]
    ans = [0.8275862068965517; 0.7586206896551725]
    return Dict("prob" => (E, F, M, γ), "ans" => ans)
end

"""
Provides the input values to replicate example 2.11 (from the book refered to 
in the introductory comments)
"""
function example_2()::Dict
    E = [1 0; 0 1]
    F = [-2; -2]
    M = [1 0; -1 0; 0 1; 0 -1]
    γ = [1; 0; 1; 0]
    ans = [1.0, 1.0]
    return Dict("prob" => (E, F, M, γ), "ans" => ans)
end

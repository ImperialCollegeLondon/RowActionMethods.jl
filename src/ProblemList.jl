
function QLR2_AN_4_6()

  problem = Model()
  @variable(nlp, x[i=1:4] >= 0, start = 0)

  @constraint(problem, - 8 +   x[1] + 2*x[2] <= 0)
  @constraint(problem, -12 + 4*x[1] +   x[2] <= 0)
  @constraint(problem, -12 + 3*x[1] + 4*x[2] <= 0)
  @constraint(problem, -8  + 2*x[3] +   x[4] <= 0)
  @constraint(problem, -8  +   x[3] + 2*x[4] <= 0)
  @constraint(problem, -5  +   x[3] +   x[4] <= 0)

  @objective(problem,  Min, 
             x[1] - x[2] - x[3] - x[1]*x[3] + x[1]*x[4] + x[2]*x[3] - x[2]*x[4])

  return problem 
end


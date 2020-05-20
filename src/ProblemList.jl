using JuMP

#All problems return their model, and an array of any variables that they 
#define. Each variables may be multivalued. 

function QLR2_AN_4_6()
  #Note; JuMP does not evaluate this QP problem to have a pos-def matrix

  problem = Model()
  #@variable(problem, x[i=1:4] >= 0, start = 0)
  @variable(problem, x[i=1:4])

  @constraint(problem, x[1] >= 0)
  @constraint(problem, x[2] >= 0)
  @constraint(problem, x[3] >= 0)
  @constraint(problem, x[4] >= 0)

  @constraint(problem, - 8 +   x[1] + 2*x[2] <= 0)
  @constraint(problem, -12 + 4*x[1] +   x[2] <= 0)
  @constraint(problem, -12 + 3*x[1] + 4*x[2] <= 0)
  @constraint(problem, -8  + 2*x[3] +   x[4] <= 0)
  @constraint(problem, -8  +   x[3] + 2*x[4] <= 0)
  @constraint(problem, -5  +   x[3] +   x[4] <= 0)

  @objective(problem,  Min, 
             x[1] - x[2] - x[3] - x[1]*x[3] + x[1]*x[4] + x[2]*x[3] - x[2]*x[4])

  return problem, x
end


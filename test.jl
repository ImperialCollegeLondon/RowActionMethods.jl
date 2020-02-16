using MathOptInterface

const MOI = MathOptInterface

model = RAM.Optimizer("Hildreth")

E = [10.0 14.0; 14.0 20.0]
F = [1.0; 2.0]

x = MOI.add_variables(model, 2)

quad = [MOI.ScalarQuadraticTerm(E[1,1],x[1],x[1]), MOI.ScalarQuadraticTerm(E[1,2],x[1],x[2]), MOI.ScalarQuadraticTerm(E[2,1],x[2],x[1]), MOI.ScalarQuadraticTerm(E[2,2],x[2],x[2])]
scal = [MOI.ScalarAffineTerm(F[i], x[i]) for i = 1:2]
obj = MOI.ScalarQuadraticFunction(scal, quad, 0.0)

MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}}(), obj)

f = MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(3.0, x[1]), MOI.ScalarAffineTerm(4.0, x[2])], 0.0)

con = MOI.add_constraint(model, f, MOI.LessThan(1.0))




using JuMP
using Ipopt 
using RowActionMethods

const RAM = RowActionMethods

model = Model(with_optimizer(() -> RAM.Optimizer("Hildreth")))
#model = Model(with_optimizer(Ipopt.Optimizer))

@variable(model, x[1:2])

E=[2 -1; -1 1]
#E=[4 0; -2 2]
F=[-1;0]
M=[-1 0; 0 -1; 3 2]
γ=[0;0;4]


@objective(model, Min, 0.5sum(x[i]*E[i,j]*x[j] for i=1:2, j=1:2) + sum(F[i]*x[i] for i=1:2))
@constraint(model, M[1,1]x[1] + M[1,2]x[2] <= γ[1])
@constraint(model, M[2,1]x[1] + M[2,2]x[2] <= γ[2])
@constraint(model, M[3,1]x[1] + M[3,2]x[2] <= γ[3])

optimize!(model)
println("x[1]: ", JuMP.value(x[1]))
println("x[2]: ", JuMP.value(x[2]))
#println("objv: ", objective_value(model))

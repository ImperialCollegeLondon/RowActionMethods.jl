[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://imperialcollegelondon.github.io/RowActionMethods.jl/dev)
![](https://github.com/ImperialCollegeLondon/RowActionMethods.jl/workflows/CI/badge.svg)
[![codecov](https://codecov.io/gh/ImperialCollegeLondon/RowActionMethods.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ImperialCollegeLondon/RowActionMethods.jl)

# Julia Row Action Optimization Solvers

This package provides a common framework for the implementation of Row Action Methods. This class of algorithms is suited to large and sparse problems. Currently only Hildreth's algorithm [[1]](#1) is implemented but other algorithms will be added in time.

The package contains a [JuMP](https://jump.dev/JuMP.jl/dev/) interface, but may also be used standalone.

<a id="1">[1]</a> Hildreth, C. (1957), A quadratic programming procedure. Naval Research Logistics, 4: 79-85. doi:10.1002/nav.3800040113

## Installation

```julia
pkg> add https://github.com/ImperialCollegeLondon/RowActionMethods.jl
```

## Usage

The following example solves a QP problem in 3 dimensions and can be formulated the same as any other JuMP problem. Note that the problem has the form `½(x'Ex)+F'x s.t. Mx≦γ`. Also note that RowActionMethods.jl defines `const RAM = RowActionMethods` as a shorthand.

### **Usage With JuMP**
```julia
using JuMP
using RowActionMethods

model = Model(with_optimizer(RAM.Optimizer, "Hildreth"))

@variable(model, x[1:3])

E=[2 3 2; 3 7 5; 2 5 4]
F=[4;7;3]
M=[-1 0 2; 0 -1 3; 3 2 4]
γ=[0;0;4]


@objective(model, Min, 0.5sum(x[i]*E[i,j]*x[j] for i=1:3, j=1:3) + sum(F[i]*x[i] for i=1:3))
@constraint(model, M[1,1]x[1] + M[1,2]x[2] + M[1,3]x[3]<= γ[1])
@constraint(model, M[2,1]x[1] + M[2,2]x[2] + M[2,3]x[3]<= γ[2])
@constraint(model, M[3,1]x[1] + M[3,2]x[2] + M[3,3]x[3]<= γ[3])

optimize!(model)
println("x[1]: ", JuMP.value(x[1]))
println("x[2]: ", JuMP.value(x[2]))
println("x[3]: ", JuMP.value(x[3]))
#Returning objective value is currently unsupported for the JuMP interface
#println("objv: ", objective_value(model))
```

### **Direct Usage**

```julia
using RowActionMethods

model = GetModel("Hildreth")

E=[2 3 2; 3 7 5; 2 5 4]
F=[4;7;3]
M=[-1 0 2; 0 -1 3; 3 2 4]
γ=[0;0;4]

SetObjective(model, E, F)

for i=1:3
    AddConstraint(model, M[i,:], γ[i])
end

optimize!(model)

println(objective_value(model))
println(GetVariables(model))
```

Solving performance is identical between JuMP and direct usage, however JuMP may be slightly slower when initially setting up the problem due a small overhead. Direct usage also allows for easier access to internal information, and new features will always be supported via the interal API before JuMP.

The internal API is defined in the documentation, and the general JuMP API is supported. If a useful feature of JuMP appears to be unsupported please raise an issue/pull request and it can be added.

## Stopping Conditions
The solver supports the ability to terminate optimisation based on an arbitrary condition. This is currently only properly supported in the direct API (but can be hacked into JuMP, proper support coming soon).

The algorithm must be used with at least one stopping condition. By default this will be an iteration limit of 32.

The built in stopping conditions are `IterationStop` and `TimeStop`. Each is instantiated with their limit (number of iterations and seconds respectively) and passed to `optimize!`.

To use a single stopping condition:
```julia
time_limit = TimeStop(30) #A 30 second time limit
optimize!(model, time_limit)
```

Use multiple stopping conditions by passing them as a vector:
```julia
time_limit = TimeStop(30) #A 30 second time limit
iteration_limit = IterationStop(15) #Stop after 15 iterations
optimize!(model, [time_limit, iteration_limit])
```

The ordering of the vector does not impact the solver. Custom stopping conditions can be defined easily. See the documentation for a guide on implementing these.

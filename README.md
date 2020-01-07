# Julia Row Action Optimization Solvers

This package contains a number of solvers defined as row-action methods. They share a common API, but each method will have specific problem features it benefits from over others. 

## API Example

A simple example using the basic Hildreth solver. All current solvers are designed for QP problems and are built with these 4 numbers, but this is not guaranteed in future. The variables E, F, m, and g refer to a QP problem in the form min ½(x'Ex)+F'x s.t. Mx≦g.

```julia
using HildrethSolver

#Select an empty model from those available in the package
model = GetModel(Hildreth())

#Create the problem representation by adding in your problem variables to the solver
buildmodel!(model, E, F, m, g)

#Repeatedly iterate the model until a stopping condition is met, by default 32 iterations
iterate_model!(m)
```

The API also offers more flexibility in stopping conditions. Note that (other than the iteration checker) each solver implements their own stopping conditions, and it is highly likely that a stopping condition for one solver will cause an error if applied to a different solver. There is currently no active check for this, so please ensure your types match. As an example of applying some limits:

```julia
model = GetModel(Hildreth())
buildmodel!(model, E, F, m, g)

#A single condition that checks the convergence values of the solver on each iteration, and terminates if it is met
convergence_condition = SC_HildrethConvergence(1e-7)

#A single condition that applies a new limit of 12 on the number of iterations
iteration_condition = SC_Iterations(12)

#Combines the conditions into a single value
#Note that this is an alias for calling:
#conditions = MultipleStopCondition([convergence_condition, iteration_condition])
conditions = get_SC(convergence_condition, iteration_condition)

iterate_model!(m, conditions)
```

## Defining Custom Stopping Conditions

In the previous example two conditions are defined (a convergence condition, and an iteration condition) and passed into the iteration function. If the builtin conditions are not sufficient, then a new one can be defined. This requires a new type that inherits from the `StoppingCondition` abstract type, and a new method for the `stopcondition` function that evaluates the new type.

```julia
#A new struct that must inherit from StoppingCondition
struct Custom_Convergence_Test <: StoppingCondition 
    #The contents are the limits to be tested against on each iteration
    limit_value ::Float64
end

#The test function, the model its solver, and your new condition type.
function stopcondition(model::HildrethModel,
                       condition::Custom_Convergence_Test
                      )::Bool
    return model.workingvars['λ'][0] < condition.limit_value
end
```

This new condition can now be used as if it was a builtin condition of the package:

```Julia 
model = GetModel(Hildreth())
buildmodel!(model, E, F, m, g)

#Create condition
new_condition = Custom_Convergence_Test(0.12)
iteration_condition = SC_Iterations(12)

#Combine conditions
conditions = get_SC(new_condition, iteration_condition)

#The condition will be tested during iteration
iterate_model!(m, conditions)
```

Please feel free to make a pull request of any useful stopping conditions so that they can be included by default.


## Known/Potential Issues
These problems should be addressed as the testing increases in scope, but please raise an issue/make a pull request with any issues you find.

- Method models may not account for the types returned by certain library functions (for example, LinearAlgebra.factorize). This will be addressed when a MOI interface is written to allow testing of a wider range of problems. 
- Current testing has only used a limited set of small problems, it is very likely that you will run into issues when using the package for slightly more complex problems. Again this will be addressed when a wider range of problems are tested.


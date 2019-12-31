# Julia Row Action Optimization Solvers

This package contains a number of solvers defined as row-action methods. They share a common API, but each method will have specific problem features it benefits from over others. 

## API Example

A simple example using the basic Hildreth solver. The variables E, F, m, and g represent the vectors/matrices within a standard QP problem. All current solvers are designed for QP problems and are built with these 4 numbers, but this is not guaranteed in future. 
```julia
using HildrethSolver

#Select an empty model from those available in the package
model = Optimizer(Hildreth())

#Create the problem representation by adding in your problem variables to the solver
buildmodel!(model, E, F, m, g)

#Repeatedly iterate the model until a stopping condition is met, by default 32 iterations
iterate_model!(m)
```

The API also offers more flexibility in stopping conditions. Note that (other than the iteration checker) each solver implements their own stopping conditions, and it is highly likely that a stopping condition for one solver will cause an error if applied to a different solver. There is currently no active check for this, so please ensure your types match. As an example of applying some limits:

```julia
model = Optimizer(Hildreth())
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

As can be seen, two new conditions are defined, combined into a single variable, and passed into the iteration function. If the builtin conditions are not sufficient, then a new one can be defined. This requires a type and a function:

```julia
#A new struct that must inherit from StoppingCondition
struct Custom_Convergence_Test <: StoppingCondition 
    #The contents are the limits to be tested against on each iteration
    limit_value ::Float64
end

#A test function, it must take the model type of the solver it should act on, and the condition struct you have just defined.
#It must return a bool. Note that the use of HildrethModel and the condition being checked is just for example. 
function Custom_Convergence_Function(model::HildrethModel,
                                     condition::Custom_Convergence_Test
                                    )::Bool
    return model.workingvars['λ'][0] < condition.limit_value
end
```


 


```@docs
GetModel(method::String)
optimize!(model::RAM.RAMProblem)
optimize!(model::RAM.RAMProblem, s::RAM.StoppingCondition)
optimize!(model::RAM.RAMProblem{T,F}, conditions::Vector{S}) where {T,F,S<:RAM.StoppingCondition}
```

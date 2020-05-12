
#TODO Add sparse vectors in Constraint entry
mutable struct ConstraintEntry{T}
    func::Vector{T}
    lim::T
end

#TODO this should be reworked to contain the problem, rather than the problem containing this

"""
    setconstraint!(model::ModelFormulation, M_row::Vector{T}, lim::T) where T

Adds a less-than constraint to the model.

Returns a UID value that the solver can use to map to the value if modifying,
viewing, or deleting the values. UID is based off the number of constraints
that have ever been added, not the current number.
"""
function addconstraint!(model::ModelFormulation, M_row::Vector{T}, lim::T)::Int where T
    !validconstraint(model, M_row, lim) && error("Invalid constraint, have you added an objective?")
    
    #Ensures unique constraint index
    new_index = model.status.max_constraint_index + 1

    push!(model.status.constraints, new_index => ConstraintEntry(M_row, lim))

    #Update largest index
    model.status.max_constraint_index = new_index
    model.status.constraint_count += 1
    return new_index
end

"""
    validconstraint(model::ModelFormulation, row::Vector{T}, lim::T)::Bool where T

Return a bool to indicate if a constraint is currently valid for the model.

Checks if the number of entries is equal to the number of registered variables.
"""
function validconstraint(model::ModelFormulation, row::Vector{T}, lim::T)::Bool where T
    return size(row)[1] == model.status.variable_count
end

"""
    extendcontraints(model::ModelFormulation)
    
Append a new value to the end of each constraint in model. Assumes that a new
variable being added is added at the end.
"""
function extendcontraints(model::ModelFormulation)
    #TODO add type parameters for appending
    for con in values(model.status.constraints)
        append!(con.func, 0.0)
    end
end

"""
    shrinkconstraints(model::ModelFormulation, index::Int)
    
Removes the entry in each constraint at index. Also removes any constraints with
no non-zero coefficients. 

Should not be called directly, as this will result in an inconsistent internal
state.
"""
function shrinkconstraints(model::ModelFormulation, index::Int)
    for con in model.status.constraints
        deleteat!(con.func, index)
    end
    filter!(c -> !check_emptyconstraint(c.second), model.status.constraints)
end

"""
    delete_variable!(model::ModelFormulation, index::Int)

Performs actions to remove constraint entries and modify objective
function to remove all references to a variable index, while maintaining
a consistent internal state. 

Note that the act of removing a variable from the objective function is
currently expensive when using a very large number of variables. This shouldn't
be an issue unless changing variables at a very regular interval.
TODO: This will likely improve when moving to sparse matrices.
"""
function delete_variable!(model::ModelFormulation, index::Int)
    model.status.variable_count -= 1
    shrinkconstraints(model, index)
    shrinkobjective(model, index)
end

"""
    check_emptyconstraint(con::ConstraintEntry)::Bool

Returns true if `con` is empty (ie, no variables are left), or if all variable
coefficients are 0.
"""
function check_emptyconstraint(con::ConstraintEntry)::Bool
    Length(con.func) == 0 && return true

    for e in con.func
        if e != 0 
            return false
        end
    end

    return true
end

"""
    get_constraintmatrix(model::ModelFormulation)::Matrix(Float64,2)

Forms the stored constraints into a transposed matrix form,
ie the value Mᵀ is returned.
"""
#TODO: Add type parameterisation
function get_constraintmatrix(model::ModelFormulation)::Array{Float64,2}
    a = hcat([j.func for j in values(model.status.constraints)]...)
    return a
end

#TODO: Add type parameterisation
function get_constraintvector(model::ModelFormulation)::Vector{Float64}
    return [j.lim for j in values(model.status.constraints)]
end

#TODO Add `reformulate` option to the model to indicate that the dual no longer represents the
#constraints
function delete_constraint!(model::ModelFormulation, con_index::Int)
    !haskey(model.status.constraints, con_index) && error("Invalid constraint identifier")
    delete!(model.status.constraints, con_index)
    model.status.constraint_count -= 1
end

function edit_constraint_coefficient!(model::ModelFormulation, con_index::Int, var_index::Int, val::Float64)
    !haskey(model.status.constraints, con_index) && error("Invalid constraint identifier")
    !(1 <= var_index <= model.status.variable_count)  && error("Invalid variable identifier")
    model.status.constraints[con_index].func[var_index] = val
end

function edit_constraint_constant!(model::ModelFormulation, con_index::Int, val::Float64)
    !haskey(model.status.constraints, con_index) && error("Invalid constraint identifier")
    model.status.constraints[con_index].lim = val
end

#TODO parametric type
function empty_model_status(model::ModelFormulation)
    s = model.status
    println(s)
    empty = s.constraints == OrderedDict{Int,ConstraintEntry{Float64}}() &&
            s.variable_count == 0 &&
            s.max_constraint_index == 0 &&
            s.constraint_count == 0 &&
            s.termination_condition == RAM_OPTIMIZE_NOT_CALLED &&
            s.iterations == 0
   return empty
end

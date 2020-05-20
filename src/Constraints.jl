export AddConstraint

#TODO this should be reworked to contain the problem, rather than the problem containing this

"""
    setconstraint!(model::ModelFormulation, M_row::Vector{T}, lim::T) where T

Adds a less-than constraint to the model.

Returns a UID value that the solver can use to map to the value if modifying,
viewing, or deleting the values. UID is based off the number of constraints
that have ever been added, not the current number.
"""
(AddConstraint(model::RAMProblem{T}, M_row::Vector, lim)::Int) where T =
    AddConstraint(model, convert(Vector{T}, M_row), convert(T, lim))

function AddConstraint(model::RAMProblem{T}, M_row::Vector{T}, lim::T)::Int where T
    !validconstraint(model, M_row, lim) && error("Invalid constraint, have you added an objective?")
    
    #Ensures unique constraint index
    new_index = model.max_constraint_index + 1

    push!(model.constraints, new_index => ConstraintEntry(M_row, lim))

    #Update largest index
    model.max_constraint_index = new_index
    model.constraint_count += 1
    return new_index
end

"""
    validconstraint(model::ModelFormulation, row::Vector{T}, lim::T)::Bool where T

Return a bool to indicate if a constraint is currently valid for the model.

Checks if the number of entries is equal to the number of registered variables.
"""
function validconstraint(model::RAMProblem, row::Vector{T}, lim::T)::Bool where T
    return size(row)[1] == model.variable_count
end

"""
    extendcontraints(model::ModelFormulation)
    
Append a new value to the end of each constraint in model. Assumes that a new
variable being added is added at the end.
"""
function extendcontraints(model::RAMProblem)
    #TODO add type parameters for appending
    for con in values(model.constraints)
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
function shrinkconstraints(model::RAMProblem, index::Int)
    for con in model.constraints
        deleteat!(con.func, index)
    end
    filter!(c -> !check_emptyconstraint(c.second), model.constraints)
end

#TODO: Move to a more appropriate place
#TODO: This will likely improve when moving to sparse matrices.
"""
    delete_variable!(model::ModelFormulation, index::Int)

Performs actions to remove constraint entries and modify objective
function to remove all references to a variable index, while maintaining
a consistent internal state. 

Note that the act of removing a variable from the objective function is
currently expensive when using a very large number of variables. This shouldn't
be an issue unless changing variables at a very regular interval.
"""
function delete_variable!(model::RAMProblem, index::Int)
    SupportsVariableDeletion(model.method) || error("Method does not support variable deletion")
    model.variable_count -= 1
    shrinkconstraints(model, index)
    DeleteVariable(model.method, index)
    DeleteSparse(model, index)
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
    GetConstraintMatrix(model::ModelFormulation)::Matrix(Float64,2)

Forms the stored constraints into a transposed matrix form,
ie the value Mᵀ is returned.
"""
(GetConstraintMatrix(model::RAMProblem{T})::Matrix{T}) where T = 
    reduce(hcat, [j.func for j in values(model.constraints)])

(GetConstraintVector(model::RAMProblem{T})::Vector{T}) where T =
    [j.lim for j in values(model.constraints)]

#TODO Add `reformulate` option to the model to indicate that the dual no longer represents the
#constraints
function delete_constraint!(model::RAMProblem, con_index::Int)
    !haskey(model.constraints, con_index) && error("Invalid constraint identifier")
    delete!(model.constraints, con_index)
    model.constraint_count -= 1
end

function edit_constraint_coefficient!(model::RAMProblem, con_index::Int, var_index::Int, val::Float64)
    !haskey(model.constraints, con_index) && error("Invalid constraint identifier")
    !(1 <= var_index <= model.variable_count)  && error("Invalid variable identifier")
    model.constraints[con_index].func[var_index] = val
end

function edit_constraint_constant!(model::RAMProblem, con_index::Int, val::Float64)
    !haskey(model.constraints, con_index) && error("Invalid constraint identifier")
    model.constraints[con_index].lim = val
end

#TODO parametric type
function is_model_empty(model::RAMProblem)
    return model.constraints == OrderedDict{Int,ConstraintEntry{Float64}}() &&
           model.variable_count == 0 &&
           model.max_constraint_index == 0 &&
           model.constraint_count == 0 &&
           model.termination_condition == RAM_OPTIMIZE_NOT_CALLED &&
           model.iterations == 0
end

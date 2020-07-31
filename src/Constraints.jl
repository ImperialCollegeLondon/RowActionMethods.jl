export AddConstraint

(AddConstraint(model::RAMProblem{T}, M_row::Vector, lim::Number)::Int) where T =
    AddConstraint(model, sparse(convert(Vector{T}, M_row)), convert(T, lim))

function AddConstraint(model::RAMProblem{T}, M_row::SparseVector{T}, lim::T)::Int where T
    !validconstraint(model, M_row, lim) && error("Invalid constraint, have you added an objective?")
    c = model.constraints
    
    if c.constraint_count == 0
        c.Functions = M_row
        c.Limits = sparse([lim])
    else
        c.Functions = hcat(c.Functions, M_row)
        c.Limits = vcat(c.Limits, lim)
    end

    c.constraint_count += 1

    #Ensures unique constraint index
    new_index = c.max_constraint_index + 1

    #Update largest index
    c.constraint_indexes[new_index] = c.constraint_count

    return new_index
end

"""
    GetConstraintMatrix(model::RAMProblem{T})::Matrix{T} where T

Return the constraint matrix.
"""
(GetConstraintMatrixTransposed(model::RAMProblem{T})::SparseMatrixCSC{T}) where T = 
    model.constraints.Functions

(GetConstraintVector(model::RAMProblem{T})::SparseVector{T}) where T = 
    model.constraints.Limits

#TODO parametric type
function is_model_empty(model::RAMProblem)
    return ((!isdefined(model.constraints, :Functions) &&
             !isdefined(model.constraints, :Limits)) || 
            (isempty(model.constraints.Functions) &&
             isempty(model.constraints.Limits))) &&
           model.variable_count == 0 &&
           model.constraints.max_constraint_index == 0 &&
           model.constraints.constraint_count == 0 &&
           model.status == OPTIMIZE_NOT_CALLED() &&
           model.iterations == 0
end

"""
    validconstraint(model::ModelFormulation, row::Vector{T}, lim::T)::Bool where T

Return a bool to indicate if a constraint is currently valid for the model.

Checks if the number of entries is equal to the number of registered variables.
"""
function validconstraint(model::RAMProblem, row::SparseVector{T}, lim::T)::Bool where T
    return size(row)[1] == model.variable_count
end

#TODO all functions below are now broken
#=

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


#TODO this needs thorough testing, not at all confident in it now
function delete_constraint!(model::RAMProblem, con_index::Int)
    RAM.SupportsDeleteConstraint(model) || error("Solver does not support constraint deletion")
    !haskey(model.constraints, con_index) && error("Invalid constraint identifier")
    delete!(model.constraints, con_index)
    removed = model.constraint_index[con_index]
    delete!(model.constraint_indexes, con_index)
    map(x=>y -> y > removed ? x=>y-1 : x=>y, model.constraint_indexes)
    model.constraint_count -= 1
    RAM.DeleteConstraint(model, removed)
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
=#

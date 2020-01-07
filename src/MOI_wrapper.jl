import MathOptInterface

const MOI = MathOptInterface
const RAM = RowActionMethods

mutable struct Optimizer <: MOI.AbstractOptimizer

    method::RAM.RowActionMethod
    inner_model::RAM.ModelFormulation

    function Optimizer(method::String;kwargs...)
        model = new()

        if !haskey(RAM.method_mapping, method) 
            throw(ArgumentError("Invalid row action method specified, valid" *
                                " methods are: $(keys(RAM.method_mapping))"))
        end

        model.method = RAM.method_mapping[method]

        for (key, val) in kwargs
            MOI.set(model, MOI.RawParameter(String(key)), val)
        end

        #MOI.empty!(model)

        return model
    end
end

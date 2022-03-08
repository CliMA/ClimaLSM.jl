module ClimaLSM
using UnPack
using DocStringExtensions
using ClimaCore
import ClimaCore: Fields
include("SharedUtilities/Domains.jl")
using .Domains
include("SharedUtilities/models.jl")
export  make_interactions_update_aux


"""
     AbstractLandModel{FT} <: AbstractModel{FT} 

An abstract type for all land model types, which are used
to simulated multiple land surface components as
a single system. Standalone component runs do not require
this interface and it should not be used for that purpose.

Many methods taking an argument of type `AbstractLandModel` are
extensions of functions defined for `AbstractModel`s. 
There are default methods that apply for all `AbstractLandModel`s,
including `make_update_aux`, `make_ode_function, `make_rhs`, 
`initialize_prognostic`, `initialize_auxiliary`, `initialize`,
and `coordinates`.

Methods which dispatch on a specific type of AbstractLandModel
include any function involving interactions between components,
as these interactions depend on the components in the land model 
and the versions of these component models being used.
"""
abstract type AbstractLandModel{FT} <: AbstractModel{FT} end


function initialize(land::AbstractLandModel)
    components = land_components(land)
    coords_list = map(components) do (component)
        Domains.coordinates(getproperty(land, component))
    end
    Y_state_list = map(zip(components, coords_list)) do (component, coords)
        getproperty(
            initialize_prognostic(getproperty(land, component), coords),
            component,
        )
    end
    p_state_list = map(zip(components, coords_list)) do (component, coords)
        getproperty(
            initialize_auxiliary(getproperty(land, component), coords),
            component,
        )
    end
    p_interactions = initialize_interactions(land)

    Y = ClimaCore.Fields.FieldVector(;
        NamedTuple(zip(components, Y_state_list))...,
    )
    p = ClimaCore.Fields.FieldVector(;
        p_interactions...,
        NamedTuple(zip(components, p_state_list))...,
    )
    coords = ClimaCore.Fields.FieldVector(;
        NamedTuple(zip(components, coords_list))...,
    )
    return Y, p, coords
end

function make_update_aux(land::AbstractLandModel)
    interactions_update_aux! = make_interactions_update_aux(land)
    components = land_components(land)
    update_aux_function_list =
        map(x -> make_update_aux(getproperty(land, x)), components)
    function update_aux!(p, Y, t)
        for f! in update_aux_function_list
            f!(p, Y, t)
        end
        interactions_update_aux!(p, Y, t) # doing this second because it needs p.soil.ψ!!
        # What if some p need other p??
        
    end
    return update_aux!
end

function make_ode_function(land::AbstractLandModel)
    components = land_components(land)
    rhs_function_list = map(x -> make_rhs(getproperty(land, x)), components)
    update_aux! = make_update_aux(land)
    function ode_function!(dY, Y, p, t)
        update_aux!(p, Y, t)
        for f! in rhs_function_list
            f!(dY, Y, p, t)
        end
    end
    return ode_function!
end

### Concrete types of AbstractLandModels
### and associated methods
include("./root_soil_model.jl")

end

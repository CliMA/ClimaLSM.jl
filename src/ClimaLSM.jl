module ClimaLSM
using UnPack
using DocStringExtensions
using ClimaCore
import ClimaCore: Fields
include("SharedUtilities/Domains.jl")
using .Domains
include("SharedUtilities/models.jl")
include("Carbon/Carbon.jl")
using .Carbon

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

ClimaLSM.name(::AbstractLandModel) = :land

function initialize(land::AbstractLandModel{FT}) where {FT}
    components = land_components(land)
    coords_list = map(components) do (component)
        Domains.coordinates(getproperty(land, component))
    end
    Y_state_list = map(zip(components, coords_list)) do (component, coords)
        zero_state = map(_ -> zero(FT), coords)
        getproperty(
            initialize_prognostic(getproperty(land, component), zero_state),
            component,
        )
    end
    p_state_list = map(zip(components, coords_list)) do (component, coords)
        zero_state = map(_ -> zero(FT), coords)
        getproperty(
            initialize_auxiliary(getproperty(land, component), zero_state),
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
        interactions_update_aux!(p, Y, t) # this has to come last if it uses p.component.value!!
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

"""
    land_components(land::AbstractLandModel)

Returns the component names of the `land` model, by calling
`propertynames(land)`.
"""
land_components(land::AbstractLandModel) = propertynames(land)

"""
    initialize_interactions(land::AbstractLandModel) end

Initializes interaction variables, which are a type of auxiliary
variable, to empty objects of the correct type
for the model. 

This function should be called during `initialize_auxiliary`,
and the method will change depending on the type of land model,
and potentially due to the type of the component models. 

This is a stub which is extended to concrete LSM models.

"""
function initialize_interactions(land::AbstractLandModel) end

"""
    make_interactions_update_aux(land::AbstractLandModel) end

Makes and returns a function which updates the interaction variables, 
which are a type of auxiliary variable.

The `update_aux!` function returned is evaluted during the right hand
side evaluation.

This is a stub which concrete types of LSMs extend.
"""
function make_interactions_update_aux(land::AbstractLandModel) end

# Methods extended by the LSM models we support
include("SurfaceWater/Pond.jl")
using .Pond
import .Pond: surface_runoff
include("Soil/Soil.jl")
using .Soil
import .Soil: source!, boundary_fluxes
include("Vegetation/Roots.jl")
using .Roots
import .Roots: flow_out_roots

### Concrete types of AbstractLandModels
### and associated methods
include("./root_soil_model.jl")
include("./pond_soil_model.jl")
end

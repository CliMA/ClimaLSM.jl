include("Soil/Soil.jl")
using .Soil
import .Soil: source
include("Vegetation/Roots.jl")
using .Roots
import .Roots: ground_area_flux_out_roots

export RootSoilModel
"""
    struct RootSoilModel{
        FT,
        SM <: Soil.AbstractSoilModel{FT},
        VM <: Roots.AbstractVegetationModel{FT},
    } <: AbstractLandModel{FT}
        soil::SM
        vegetation::VM
    end

A concrete type of land model used for simulating systems with a 
vegetation and a soil component.
$(DocStringExtensions.FIELDS)
"""
struct RootSoilModel{
    FT,
    SM <: Soil.AbstractSoilModel{FT},
    VM <: Roots.AbstractVegetationModel{FT},
} <: AbstractLandModel{FT}
    "The soil model to be used"
    soil::SM
    "The vegetation model to be used"
    vegetation::VM
end


"""
    RootSoilModel{FT}(;
                      soil_model_type::Type{SM},
                      soil_args::NamedTuple = (;),
                      vegetation_model_type::Type{VM},
                      vegetation_args::NamedTuple = (;),
                      ) where {FT,
                               SM <: Soil.AbstractSoilModel{FT},
                               VM <: Roots.AbstractVegetationModel{FT}}

A constructor for the `RootSoilModel`, which takes in the concrete model
type and required arguments for each component, constructs those models,
and constructs the `RootSoilModel` from them.

Each component model is constructed with everything it needs to be stepped
forward in time, including boundary conditions, source terms, and interaction
terms.
"""
function RootSoilModel{FT}(;
    soil_model_type::Type{SM},
    soil_args::NamedTuple = (;),
    vegetation_model_type::Type{VM},
    vegetation_args::NamedTuple = (;),
) where {
    FT,
    SM <: Soil.AbstractSoilModel{FT},
    VM <: Roots.AbstractVegetationModel{FT},
}

    #These may be passed in, or set, depending on use scenario
    boundary_fluxes = FluxBC{FT}(FT(0.0), FT(0.0))
    transpiration = PrescribedTranspiration{FT}((t::FT) -> FT(0.0))

    ##These should always be set by the constructor.
    sources = (RootExtraction{FT}(),)
    root_extraction = PrognosticSoilPressure{FT}()

    soil = soil_model_type(;
        boundary_conditions = boundary_fluxes,
        sources = sources,
        soil_args...,
    )
    vegetation = vegetation_model_type(;
        root_extraction = root_extraction,
        transpiration = transpiration,
        vegetation_args...,
    )
    args = (soil, vegetation)
    return RootSoilModel{FT, typeof.(args)...}(args...)
end

"""
    land_components(land::RootSoilModel)

Returns the component names of the `RootSoilModel`.
"""
land_components(land::RootSoilModel) = (:soil, :vegetation)

"""
    initialize_interactions(
        land::RootSoilModel{FT, SM, RM},
    ) where {FT, SM <: Soil.RichardsModel{FT}, RM <: Roots.RootsModel{FT}}

Initializes interaction variables, which are a type of auxiliary
variable, to empty objects of the correct type
for the model. 

This function should be called during `initialize_auxiliary`,
and the method will change depending on the type of land model,
and potentially due to the type of the component models. 
"""
function initialize_interactions(#Do we want defaults, for land::AbstractLandModel?
    land::RootSoilModel{FT, SM, RM},
) where {FT, SM <: Soil.RichardsModel{FT}, RM <: Roots.RootsModel{FT}}

    soil_coords = land.soil.coordinates
    return (root_extraction_source = similar(soil_coords), net_root_extraction_flux = zero(FT))
end

"""
    make_interactions_update_aux(
        land::RootSoilModel{FT, SM, RM},
    ) where {FT, SM <: Soil.RichardsModel{FT}, RM <: Roots.RootsModel{FT}}

Makes and returns a function which updates the interaction variables, 
which are a type of auxiliary variable.

The `update_aux!` function returned is evaluted during the right hand
side evaluation.
"""
function make_interactions_update_aux(#Do we want defaults, for land::AbstractLandModel?
    land::RootSoilModel{FT, SM, RM},
) where {FT, SM <: Soil.RichardsModel{FT}, RM <: Roots.RootsModel{FT}}
    function update_aux!(p, Y, t)
        @unpack a_root, b_root, K_max_root, =
            land.vegetation.param_set
        gf = ground_area_flux.(
            land.soil.coordinates,# we have to use coordinates here rather than root depth array to avoid broadcast error? 
            land.vegetation.domain.compartment_heights[1],
             p.soil.ψ.* FT(9800), # Pa: ρg * meters
            Roots.theta_to_p(Y.vegetation.theta[1]),
            a_root,
            b_root,
            K_max_root,
        ) .* land.vegetation.param_set.root_distribution_function.(land.soil.coordinates) # here too
        @. p.root_extraction_source =  -gf * land.vegetation.param_set.RAI
        p.net_root_extraction_flux = sum(gf)
    end
    return update_aux!
end


## Extending methods of the Roots and Soil Model
## TBD if these should be linked more explicitly.
"""
   PrognosticSoilPressure{FT} <: Roots.AbstractRootExtraction{FT}

The concrete type of root extraction model, used for dispatch when computing
the right hand side of the vegetation model when soil is also modeled
prognostically.

This is paired with the source term `RootExtraction`, which returns 
the flow of water between roots and soil in units of 1/sec, 
rather than moles/sec, as needed by the soil model.
"""
struct PrognosticSoilPressure{FT} <: Roots.AbstractRootExtraction{FT} end

"""
    Roots.flow_out_roots(
        re::PrognosticSoilPressure{FT},
        model::Roots.RootsModel{FT},
        Y::ClimaCore.Fields.FieldVector,
        p::ClimaCore.Fields.FieldVector,
        t::FT,
    )::FT where {FT}

An extension of the function which returns the
net flow of water between the
roots and the soil, when both soil and plant
hydraulics are modeled prognostically.

It is computed by summing the flow of water between
roots and soil at each soil layer.
"""
function Roots.ground_area_flux_out_roots(
    re::PrognosticSoilPressure{FT},
    model::Roots.RootsModel{FT},
    Y::ClimaCore.Fields.FieldVector,
    p::ClimaCore.Fields.FieldVector,
    t::FT,
)::FT where {FT}
        return p.net_root_extraction_flux
end

"""
    RootExtraction{FT} <: Soil.AbstractSoilSource{FT}

The concrete type of root extraction model, used for dispatch when computing
the flow of water between soil and the roots as a source term for the
soil model (1/s).

This is paired with the root model type `PrognosticSoilPressure`, which 
returns the flow of water between roots and soil in units of moles/sec, 
as needed by the root model.
"""
struct RootExtraction{FT} <: Soil.AbstractSoilSource{FT} end

"""
    Soil.source(src::RootExtraction{FT}, Y, p) where {FT}

An extension of the function which computes source terms for the 
soil model; this method returns the water loss or gain due
to roots when a plant hydraulic prognostic model is included.
"""
function Soil.source(src::RootExtraction{FT}, Y, p) where {FT}
    return p.root_extraction_source
end

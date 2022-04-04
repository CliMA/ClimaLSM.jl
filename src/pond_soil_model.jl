export LandHydrology, infiltration_capacity, infiltration_at_point

"""
    struct LandHydrology{
        FT,
        SM <: Soil.AbstractSoilModel{FT},
        SW <: Pond.AbstractSurfaceWaterModel{FT},
    } <: AbstractLandModel{FT}

A concrete type of land model used for simulating systems with a 
soil and surface water component.
$(DocStringExtensions.FIELDS)"""
struct LandHydrology{
    FT,
    SM <: Soil.AbstractSoilModel{FT},
    SW <: Pond.AbstractSurfaceWaterModel{FT},
} <: AbstractLandModel{FT}
    "The soil model"
    soil::SM
    "The surface water model"
    surface_water::SW
end

"""
    LandHydrology{FT}(;
        land_args::NamedTuple = (;),
        soil_model_type::Type{SM},
        soil_args::NamedTuple = (;),
        surface_water_model_type::Type{SW},
        surface_water_args::NamedTuple = (;),
    ) where {
        FT,
        SM <: Soil.AbstractSoilModel{FT},
        SW <: Pond.AbstractSurfaceWaterModel{FT},
    }
A constructor for the `LandHydrology` model, which takes in the concrete model
type and required arguments for each component, constructs those models,
and constructs the `LandHydrology` from them.

Each component model is constructed with everything it needs to be stepped
forward in time, including boundary conditions, source terms, and interaction
terms.

Additional arguments, like parameters and driving atmospheric data, are passed
in as `land_args`.
"""
function LandHydrology{FT}(;
    land_args::NamedTuple = (;),
    soil_model_type::Type{SM},
    soil_args::NamedTuple = (;),
    surface_water_model_type::Type{SW},
    surface_water_args::NamedTuple = (;),
) where {
    FT,
    SM <: Soil.AbstractSoilModel{FT},
    SW <: Pond.AbstractSurfaceWaterModel{FT},
}

    @unpack precip = land_args

    sources = ()
    surface_runoff = PrognosticRunoff{FT}(precip)
    boundary_conditions = RunoffBC{FT}()

    soil = soil_model_type(;
        boundary_conditions = RunoffBC{FT}(),
        sources = sources,
        soil_args...,
    )
    surface_water = surface_water_model_type(;
        runoff = surface_runoff,
        surface_water_args...,
    )
    args = (soil, surface_water)
    return LandHydrology{FT, typeof.(args)...}(args...)
end

"""
    function initialize_interactions(
        land::LandHydrology{FT, SM, SW},
    ) where {FT, SM <: Soil.RichardsModel{FT}, SW <: Pond.PondModel{FT}}

Initializes the interaction auxiliary variable `p.soil_infiltration`.
"""
function initialize_interactions(
    land::LandHydrology{FT, SM, SW},
) where {FT, SM <: Soil.RichardsModel{FT}, SW <: Pond.PondModel{FT}}
    # Do we do coordinats(surface) and compute it another time (new instance of space)
    # or store it? Then we get (x,y)
    # Do we do similar(soil coordinates) and not have to store or compute anything again,
    # but then get (x,y,z = z of some center layer)?
    # Which is better for broadcasting later on?

    # If the surface domain is a hybrid box, is that just easier?
    # Also, if the surface has topography, how else would we carry along z info?
    # it would only live in the soil model as is
    # That may also be more consistent with what we do with Point{FT}(z_sfc)
    surface_coords = coordinates(land.surface_water)#In non-column case, would this be level(soil_coords)?
    return (soil_infiltration = similar(surface_coords),)
end

#=
If there is a pond present, flux BC should be -i_c

If there is no pond present (we won't model the standing surface water in full LSM)

 _______________________________P-E______________________________________________
                      |   P- E >0 (into soil)      |   P-E <0 (out of soil)
 |________________________________________________________________________________
 | Into Soil, i_c >0  |   min(P-E, i_c)            |   P-E -> min(P-E, i_c) works
I|________________________________________________________________________________
 |Out of soil, i_c <0 |   i_c; min(P-E, i_c) works |   Should be E? rare? unclear?
 |________________________________________________________________________________
=#

"""
    infiltration_at_point(η::FT, i_c::FT, P::FT)   

Returns the infiltration given pond height η, infiltration capacity,
and precipitation. 


This is defined such that positive means into soil.
"""
infiltration_at_point(η::FT, i_c::FT, P::FT) where {FT} =
    η > eps(FT) ? i_c : min(i_c, P)

"""
    function infiltration_capacity(
        model::Soil.AbstractSoilModel{FT},
        Y::ClimaCore.Fields.FieldVector,
        p::ClimaCore.Fields.FieldVector,
    ) where {FT}

Function which computes the infiltration capacity of the soil based on
soil characteristics, moisture levels, and pond height.

Defined such that positive means into soil.
"""
function infiltration_capacity(
    model::Soil.AbstractSoilModel{FT},
    Y::ClimaCore.Fields.FieldVector,
    p::ClimaCore.Fields.FieldVector,
) where {FT}
    N = length(axes(p.ψ).center_local_geometry)
    Δz = last(axes(p.ψ).face_local_geometry.WJ)
    z_surf = Δz + last(parent(model.coordinates))
    ψ, K =
        ClimaCore.Operators.getidx.(
            [p.ψ, p.K],
            Ref(ClimaCore.Operators.Interior()),
            N,
        )
        #approx -k∇(psi +z) by setting face value and looking up center value
        #K(top center)* ( psi(face) - psi(center))/dz
    return @. (K * (max(0.0, Y.surface_water.η) - ψ + Δz) / Δz)
end

"""
    function make_interactions_update_aux(
        land::LandHydrology{FT, SM, SW},
    ) where {FT, SM <: Soil.RichardsModel{FT}, SW <: Pond.PondModel{FT}}

A method which makes a function; the returned function 
updates the auxiliary variable `p.soil_infiltration`, which
is needed for both the boundary condition for the soil model and the source
term (runoff) for the surface water model.

This function is called each ode function evaluation.
"""
function make_interactions_update_aux(
    land::LandHydrology{FT, SM, SW},
) where {FT, SM <: Soil.RichardsModel{FT}, SW <: Pond.PondModel{FT}}
    function update_aux!(p, Y, t)
        i_c = infiltration_capacity(land.soil, Y, p.soil)
        @. p.soil_infiltration =
            -infiltration_at_point(
                Y.surface_water.η,
                i_c,
                -land.surface_water.runoff.precip(t),
            )
    end
    return update_aux!
end

"""
    PrognosticRunoff{FT} <: Pond.AbstractSurfaceRunoff{FT}

Concrete type of `Pond.AbstractSurfaceRunoff` for use in LSM models,
where precipitation is passed in, but infiltration is computed
prognostically.

This is paired with `Soil.RunoffBC`: both are used at the same time,
ensuring the infiltration used for the boundary condition of soil
is also used to compute the runoff for the surface water.

"""
struct PrognosticRunoff{FT} <: Pond.AbstractSurfaceRunoff{FT}
    precip::Function
end

"""
    function Pond.surface_runoff(
        runoff::PrognosticRunoff{FT},
        Y::ClimaCore.Fields.FieldVector,
        p::ClimaCore.Fields.FieldVector,
        t::FT,
    ) where {FT}

Extension of the `Pond.surface_runoff` function, which computes the surface runoff, for use in an LSM when the runoff is determined
prognostically.
"""
function Pond.surface_runoff(
    runoff::PrognosticRunoff{FT},
    Y::ClimaCore.Fields.FieldVector,
    p::ClimaCore.Fields.FieldVector,
    t::FT,
) where {FT}
    return @. -(runoff.precip(t) - p.soil_infiltration)
end

"""
    RunoffBC{FT} <: Soil.AbstractSoilBoundaryConditions{FT}

Concrete type of `Soil.AbstractSoilBoundaryConditions` for use in LSM models,
where precipitation is passed in, but infiltration is computed
prognostically. This infiltration is then used to set an upper
boundary condition for the soil.


This is paired with `Pond.PrognosticRunoff`: both are used at the same 
time,
ensuring that the infiltration used for the boundary condition of soil
is also used to compute the runoff for the surface water.
"""
struct RunoffBC{FT} <: Soil.AbstractSoilBoundaryConditions{FT} end

"""
    function Soil.boundary_fluxes(
        bc::RunoffBC{FT},
        p::ClimaCore.Fields.FieldVector,
        t::FT,
    ) where {FT}

Extension of the `Soil.boundary_fluxes` function, which
 returns the water volume boundary flux for the soil, by
returning the soil infiltration (computed each step and stored in
`p.soil_infiltration`) as a top flux, and 
assuming no flux at the bottom of the soil. This is for use in an LSM.
"""
function Soil.boundary_fluxes(
    bc::RunoffBC{FT},
    p::ClimaCore.Fields.FieldVector,
    t::FT,
) where {FT}
    return p.soil_infiltration[1], FT(0.0)
end

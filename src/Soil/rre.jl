"""
    RichardsParameters{FT <: AbstractFloat}

A struct for storing parameters of the `RichardModel`.
$(DocStringExtensions.FIELDS)
"""
struct RichardsParameters{FT <: AbstractFloat}
    "The porosity of the soil (m^3/m^3)"
    ν::FT
    "The van Genuchten parameter α (1/m)"
    vg_α::FT
    "The van Genuchten parameter n"
    vg_n::FT
    "The van Genuchten parameter m"
    vg_m::FT
    "The saturated hydraulic conductivity (m/s)"
    Ksat::FT
    "The specific storativity (1/m)"
    S_s::FT
    "The residual water fraction (m^3/m^3"
    θ_r::FT
end

"""
    RichardsModel

A model for simulating the flow of water in a porous medium
by solving the Richardson-Richards Equation.

$(DocStringExtensions.FIELDS)
"""
struct RichardsModel{FT, PS, D, C, BC, S} <: AbstractSoilModel{FT}
    "the parameter set"
    param_set::PS
    "the soil domain, using ClimaCore.Domains"
    domain::D
    "the domain coordinates"
    coordinates::C
    "the boundary conditions, of type AbstractSoilBoundaryConditions"
    boundary_conditions::BC
    "A tuple of sources, each of type AbstractSoilSource"
    sources::S
end

"""
    RichardsModel{FT}(;
        param_set::RichardsParameters{FT},
        domain::D,
        boundary_conditions::AbstractSoilBoundaryConditions{FT},
        sources::Tuple,
    ) where {FT, D}

A constructor for a `RichardsModel`.
"""
function RichardsModel{FT}(;
    param_set::RichardsParameters{FT},
    domain::D,
    boundary_conditions::AbstractSoilBoundaryConditions{FT},
    sources::Tuple,
) where {FT, D}
    coords = coordinates(domain)
    args = (param_set, domain, coords, boundary_conditions, sources)
    RichardsModel{FT, typeof.(args)...}(args...)
end


"""
    make_rhs(model::RichardsModel)

An extension of the function `make_rhs`, for the Richardson-
Richards equation. 

This function creates and returns a function which computes the entire
right hand side of the PDE for `ϑ_l`, and updates `dY.soil.ϑ_l` in place
with that value.

This has been written so as to work with Differential Equations.jl.
"""
function ClimaLSM.make_rhs(model::RichardsModel)
    function rhs!(dY, Y, p, t)
        @unpack ν, vg_α, vg_n, vg_m, Ksat, S_s, θ_r = model.param_set
        top_flux_bc, bot_flux_bc =
            boundary_fluxes(model.boundary_conditions, p, t)
        z = model.coordinates.z
        interpc2f = Operators.InterpolateC2F()
        gradc2f_water = Operators.GradientC2F()
        divf2c_water = Operators.DivergenceF2C(
            top = Operators.SetValue(Geometry.WVector(top_flux_bc)),
            bottom = Operators.SetValue(Geometry.WVector(bot_flux_bc)),
        )
        @. dY.soil.ϑ_l =
            -(divf2c_water(-interpc2f(p.soil.K) * gradc2f_water(p.soil.ψ + z)))
        # Horizontal contributions
        horizontal_components!(dY, model.domain, model, p)

        # Source terms
        for src in model.sources
            source!(dY, src, Y, p)
        end

        # This has to come last
        dss!(dY, model.domain)
    end
    return rhs!
end


"""
   horizontal_components!(dY::ClimaCore.Fields.FieldVector,
                          domain::HybridBox,
                          model::RichardsModel,
                          p::ClimaCore.Fields.FieldVector)
Updates dY in place by adding in the tendency terms resulting from
horizontal derivative operators.

In the case of a hybrid box domain, the horizontal contributions are
computed using the WeakDivergence and Gradient operators.
"""
function horizontal_components!(
    dY::ClimaCore.Fields.FieldVector,
    domain::HybridBox,
    model::RichardsModel,
    p::ClimaCore.Fields.FieldVector,
)
    hdiv = Operators.WeakDivergence()
    hgrad = Operators.Gradient()
    @. dY.soil.ϑ_l += -hdiv(-p.soil.K * hgrad(p.soil.ψ + model.coordinates.z))
end

"""
    prognostic_vars(soil::RichardsModel)

A function which returns the names of the prognostic variables
of `RichardsModel`.
"""
ClimaLSM.prognostic_vars(soil::RichardsModel) = (:ϑ_l,)

"""
    auxiliary_vars(soil::RichardsModel)

A function which returns the names of the auxiliary variables 
of `RichardsModel`.

Note that auxiliary variables are not needed for such a simple model.
We could instead compute the conductivity and matric potential within the
rhs function explicitly, rather than compute and store them in the 
auxiliary vector `p`. We did so in this case as a demonstration.
"""
ClimaLSM.auxiliary_vars(soil::RichardsModel) = (:K, :ψ)

"""
    make_update_aux(model::RichardsModel)

An extension of the function `make_update_aux`, for the Richardson-
Richards equation. 

This function creates and returns a function which updates the auxiliary
variables `p.soil.variable` in place.

This has been written so as to work with Differential Equations.jl.
"""
function ClimaLSM.make_update_aux(model::RichardsModel)
    function update_aux!(p, Y, t)
        @unpack ν, vg_α, vg_n, vg_m, Ksat, S_s, θ_r = model.param_set
        @. p.soil.K = hydraulic_conductivity(
            Ksat,
            vg_m,
            effective_saturation(ν, Y.soil.ϑ_l, θ_r),
        )
        @. p.soil.ψ = pressure_head(vg_α, vg_n, vg_m, θ_r, Y.soil.ϑ_l, ν, S_s)
    end
    return update_aux!
end

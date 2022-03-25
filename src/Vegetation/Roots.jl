# To do: convert all units to SI, move constants out of code and use ClimaParameters,
# convert from a single plant to a bulk plant (units change, required input changes).
# we should also change the name to Vegetation or biophysics as appropriate.
module Roots
#=
    Roots

This module contains everything needed to run a vegetation model
in standalone mode.

The vegetation model is assumed to have a set of prognostic `Y` and
auxiliary `p` variables, which describe the state of the
vegetation system. The system is evolved in time by solving 
equations of the form

```

\frac{d Y}{d t} = f(Y, t, p(Y, t; \ldots);\ldots),

```

i.e. ordinary differential equations depending on state `Y`,
auxiliary functions of the state `p`, and other parameters 
represented by the ellipses. For example, `p` may represent
the transpiration rate, which must be computed each time step
based on the vegetation prognostic state, functions of time, 
like the atmosphere state, and other parameters.

Currently, only a simple plant hydraulics model is supported,
but our plan is to include much more complex representations
of the vegetation.

Addition of additional versions of vegetation
models requires defining a model type (of super type 
`AbstractVegetationModel`), and extending the methods
imported by Models.jl, as needed, for computing the
right hand side functions of the ordinary differential equations
and the functions which update auxiliary variables whenever
the right hand side is evaluated.

This code base assumes that DifferentialEquations.jl
will be used for evolving the system in time,
and that the array-like objected being stepped forward
is a `ClimaCore.Fields.FieldVector`. 

While a simple array may be sufficient for vegetation models,
the `FieldVector` type is used in order to make use of 
`ClimaCore` functionality when solving PDEs (required for
other components of Land Surface Models) and for ease of
handling multi-column models.

To simulated land surfaces with multiple components (vegetation,
soil, rivers, etc), the ClimaLSM.jl package should be used.
That package will use the methods of this function for advancing
the system forward in time, extending methods as needed to account
for interactions between components.
=#
using ClimaLSM
using ClimaCore
using UnPack
using DocStringExtensions
import ClimaCore: Fields

using ClimaLSM.Domains: AbstractVegetationDomain, RootDomain
import ClimaLSM:
    AbstractModel,
    initialize_prognostic,
    make_update_aux,
    make_rhs,
    make_ode_function,
    prognostic_vars,
    auxiliary_vars,
    initialize,
    initialize_auxiliary,
    name
export RootsModel,
    AbstractVegetationModel,
    ground_area_flux,
    ground_area_flux_out_roots,
    θ_to_p,
    p_to_θ,
    RootsParameters,
    PrescribedSoilPressure,
    PrescribedTranspiration,
    AbstractRootExtraction

"""
    AbstractVegetationModel{FT} <: AbstractModel{FT}

An abstract type for vegetation models.

Concrete types include a plant hydraulics model, but future types will
include multi-layer canopy models and possibly a big leaf model.
"""
abstract type AbstractVegetationModel{FT} <: AbstractModel{FT} end

ClimaLSM.name(::AbstractVegetationModel) = :vegetation

"""
    AbstractRootExtraction{FT <: AbstractFloat}

An abstract type for types representing different models of
water exchange between soil and plants.

Currently, only a prescribed soil pressure is supported for
standalone plant hydraulics. Use within an LSM requires
types defined within ClimaLSM, and include a prognostic
soil pressure for models with both soil and roots.

"""
abstract type AbstractRootExtraction{FT <: AbstractFloat} end

"""
    AbstractTranspiration{FT <: AbstractFloat}

An abstract type for types representing different models of
transpiration.

Currently, only a PrescribedTranspiration is supported.
"""
abstract type AbstractTranspiration{FT <: AbstractFloat} end


"""
    RootsParameters{FT <: AbstractFloat}


A struct for holding parameters of the Root Model. Eventually to be used with ClimaParameters.
$(DocStringExtensions.FIELDS)
"""
struct RootsParameters{FT <: AbstractFloat}
    "controls the shape and steepness of conductivity vs. pressure curve, for roots: unitless"
    a_root::FT
    "controls the steepness of the relative conductivity vs. pressure curve, for roots: inverse Pa"
    b_root::FT
    "controls the shape and steepness of relative conductivity vs. pressure curve, for stems: unitless"
    a_stem::FT
    "controls the steepness of the conductivity vs. pressure curve, for stems: inverse Pa"
    b_stem::FT
    "height of stem (m)"
    h_stem::FT
    #"thickness of the leaves (m)"
    #h_leaf::FT
    "maximum water conductivity in roots (m^3*m/s/Pa/m^2 conducting area)"
    K_max_root::FT
    "maximum water conductivity in stems (m^3*m/s/Pa/m^2 conducting area)"
    K_max_stem::FT
    "Leaf area index: surface area of leaves/area of ground"
    #LAI::FT
    #"root area index: cross section of roots/area of ground"
    #RAI::FT
    #"Stem area index: cross section of stem/area of ground"
    #SAI::FT
    #"Root distribution function P(z)"
    #root_distribution_function::Function
end

"""
    RootsModel{FT, PS, D, RE, T, B} <: AbstractVegetationModel{FT}

Defines, and constructs instances of, the RootsModel type, which is used
for simulation ground_area_flux of water to/from soil, along roots of different depths,
along a stem, to a leaf, and ultimately being lost from the system by
transpiration. 

This model can be used in standalone mode by prescribing the transpiration rate
and soil pressure at the root tips, or with a dynamic soil model using `ClimaLSM`.

$(DocStringExtensions.FIELDS)
"""
struct RootsModel{FT, PS, D, RE, T} <: AbstractVegetationModel{FT}
    "Parameters required by the root model"
    param_set::PS
    "The root model domain, of type `AbstractVegetationDomain`"
    domain::D
    "The root extraction model, of type `AbstractRootExtraction`"
    root_extraction::RE
    "The transpiration model, of type `AbstractTranspiration`"
    transpiration::T
end

function RootsModel{FT}(;
    param_set,
    domain::AbstractVegetationDomain{FT},
    root_extraction::AbstractRootExtraction{FT},
    transpiration::AbstractTranspiration{FT},
) where {FT}
    args = (param_set, domain, root_extraction, transpiration)
    return RootsModel{FT, typeof.(args)...}(args...)
end

"""
    function modified_logisitic_vulnerability_curve(
        a::FT,
        b::FT,
        p::FT,
    )::FT where {FT}

Computes the conductivity given parameters a, b and pressure.
"""
function modified_logisitic_vulnerability_curve(
    a::FT,
    b::FT,
    p::FT,
)::FT where {FT}
    v_c = (a + 1)/a * (1 - 1/(1+a*exp(b*p)))  
    return v_c
end

"""
    function modified_logistic_int_k_dp
        a::FT,
        b::FT,
        p::FT,
    )::FT where {FT}

Computes the integral of the modified logisitc vulnerability curve.
"""
function modified_logistic_int_k_dp(
    a::FT,
    b::FT,
    p1::FT,
    p2::FT,
)::FT where {FT}
    int_k_dp1 = (a + 1)/a * (log(a*exp(b*p1)+1))/b  
    int_k_dp2 = (a + 1)/a * (log(a*exp(b*p2)+1))/b 
    int_k_dp = int_k_dp2 - int_k_dp1
    return int_k_dp
end

"""
    prognostic_vars(model::RootsModel)

A function which returns the names of the prognostic 
variables of the `RootsModel`.
"""
prognostic_vars(model::RootsModel) = (:θ,)

"""
    function ground_area_flux(
        z1::FT,
        z2::FT,
        p1::FT,
        p2::FT,
        a::FT,
        b::FT,
        Kmax::FT
    ) where {FT}

Computes the ground_area_flux of water (volume of water/ground area/second)  given the height and pressures
at two points. Here, `a`, `b, and `Kmax` are parameters
which parameterize the hydraulic conductance of the pathway along which
the ground_area_flux occurs. `AI` is the area index relating conducting area to ground area.
"""
function ground_area_flux(
    z1::FT,
    z2::FT,
    p1::FT,
    p2::FT,
    a::FT,
    b::FT,
    K_max_root::FT
)::FT where {FT}
    ρg = FT(0.0098) # Pa/m 
    cond_area_flux = -K_max_root*modified_logistic_int_k_dp(a,b,p1,p2) - ρg*(vulnerability_curve(a,b,p1) + vulnerability_curve(a,b,p2))/2
    return cond_area_flux
end


"""
    θ_to_p(θ::FT) where {FT}

Computes the volumetric water content given pressure (p).
"""
function θ_to_p(θ::FT) where {FT}
    θ = min(θ, FT(1.0))
    θ = max(eps(FT), θ)
    p = (θ-1)*5  
    return p
end


"""
    p_to_θ(p::FT) where {FT}

Computes the pressure (p)  given the volumetric water content (θ).
"""
function p_to_θ(p::FT) where {FT}
    θ = p/5+1
    return θ
end


"""
    make_rhs(model::RootsModel)

A function which creates the rhs! function for the RootsModel.

The rhs! function must comply with a rhs function of OrdinaryDiffEq.jl.
"""
function make_rhs(model::RootsModel{FT}) where {FT}
    function rhs!(dY, Y, p, t)
        @unpack a_root,
        b_root,
        K_max_root,
        a_stem,
        b_stem,
        K_max_stem,
        h_stem = model.param_set

        z_stem, z_leaf = model.domain.compartment_heights

        p_stem = θ_to_p(Y.vegetation.θ[1])
        p_leaf = θ_to_p(Y.vegetation.θ[2])

        # Includes RAI factor
        ground_area_flux_in_stem = ground_area_flux_out_roots(model.root_extraction, model, Y, p, t)

        # Includes SAI factor
        ground_area_flux_out_stem = ground_area_flux(
            z_stem,
            z_leaf,
            p_stem,
            p_leaf,
            a_stem,
            b_stem,
            K_max_stem,
        )

        dY.vegetation.θ[1] = ground_area_flux_in_stem - ground_area_flux_out_stem
        dY.vegetation.θ[2] = ground_area_flux_out_stem - ground_area_transpiration(model, model.transpiration, t)
    end
    return rhs!
end

"""
    PrescribedSoilPressure{FT} <: AbstractRootExtraction{FT}

A concrete type used for dispatch when computing the `ground_area_flux_out_roots`,
in the case where the soil pressure at each root layer is prescribed.
"""
struct PrescribedSoilPressure{FT} <: AbstractRootExtraction{FT}
    p_soil::Function
end

"""
    PrescribedTranspiration{FT} <: AbstractTranspiration{FT}

A concrete type used for dispatch when computing the transpiration
from the leaves, in the case where transpiration is prescribed.
"""
struct PrescribedTranspiration{FT} <: AbstractTranspiration{FT}
    T::Function
end

"""
    ground_area_flux_out_roots(
        re::PrescribedSoilPressure{FT},
        model::RootsModel{FT},
        Y::ClimaCore.Fields.FieldVector,
        p::ClimaCore.Fields.FieldVector,
        t::FT,
    )::FT where {FT}

A method which computes the ground_area_flux between the soil and the stem, via the roots,
in the case of a standalone root model with prescribed soil pressure
at the root tips.

This assumes that the stem compartment is the first element of `Y.roots.θ`.
"""
function ground_area_flux_out_roots(
    re::PrescribedSoilPressure{FT},
    model::RootsModel{FT},
    Y::ClimaCore.Fields.FieldVector,
    p::ClimaCore.Fields.FieldVector,
    t::FT,
)::FT where {FT}
    @unpack a_root, b_root, K_max_root =
        model.param_set
    p_stem = θ_to_p(Y.vegetation.θ[1])
    return sum(
        ground_area_flux.(
            model.domain.root_depths,
            model.domain.compartment_heights[1],
            re.p_soil(t),
            p_stem,
            a_root,
            b_root,
            K_max_root
        ).* (vcat(model.domain.root_depths,[0.0])[2:end] - vcat(model.domain.root_depths,[0.0])[1:end-1])) 
        #.* model.param_set.root_distribution_function.(model.domain.root_depths)
    end

"""
    ground_area_transpiration(model::RootsModel{FT},
        transpiration::PrescribedTranspiration{FT},
        t::FT,
    )::FT where {FT}

A method which computes the transpiration in volume of water/ground area/second between the leaf
and the atmosphere,
in the case of a standalone root model with prescribed transpiration rate.
"""
function ground_area_transpiration(model::RootsModel{FT},
    transpiration::PrescribedTranspiration{FT},
    t::FT,
)::FT where {FT}
    return transpiration.T(t)
end

end

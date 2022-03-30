# To do: add in SoilHeatWaterModel, different boundary conditions, freeze thaw.
module Soil
#=
    Soil

This module contains everything needed to run a soil model
in standalone mode.

The soil model is assumed to have a set of prognostic `Y` and
auxiliary `p` variables, which describe the state of the
soil system. The system is evolved in time by solving 
equations of the form

```

\frac{\partial Y}{\partial t} = D(Y, t, p(Y, t; \ldots);\ldots),

```

i.e. partial (or ordinary) differential equations depending 
on state `Y`,
auxiliary functions of the state `p`, and other parameters 
represented by the ellipses. The operator `D` indicates a generic
nonlinear differential operator.  Not every model
requires auxilary variables, but these are useful for storing
quantities that are needed multiple times per right hand side
evaluation. An example would be the temperature, which is computed
from the prognostic variables of water content and internal
energy.

Currently, both the Richardson Richards Equation (RRE; hydrology alone)
and an integrated soil energy and hydrology model are supported.

Addition of additional versions of soil
models requires defining a model type (of super type 
`AbstractSoilModel`), and extending the methods
imported by Models.jl, as needed, for computing the
right hand side functions of the ordinary differential equations
and the functions which update auxiliary variables whenever
the right hand side is evaluated.

This code base assumes that DifferentialEquations.jl
will be used for evolving the system in time,
and that the array-like objected being stepped forward
is a `ClimaCore.Fields.FieldVector`. 

The `FieldVector` type is used in order to make use of 
`ClimaCore` functionality when solving PDEs (`ClimaCore` handles
all of the spatial discretization and operator functionality)
 and for ease of handling multi-column models.

To simulate land surfaces with multiple components (vegetation,
soil, rivers, etc), the ClimaLSM.jl package should be used.
That package will use the methods of this function for advancing
the system forward in time, extending methods as needed to account
for interactions between components.
=#

using ClimaLSM
using UnPack
using DocStringExtensions
using ClimaCore
import ClimaCore: Fields, Operators, Geometry, Spaces

import ClimaLSM.Domains: coordinates, Column, HybridBox
import ClimaLSM:
    AbstractModel,
    make_update_aux,
    make_rhs,
    prognostic_vars,
    auxiliary_vars,
    name
export RichardsModel,
    SoilEnergyHydrology,
    RichardsParameters,
    HeatParameters,
    boundary_fluxes,
    FluxBC,
    PrecipFreeDrainage,
    RootExtraction,
    AbstractSoilModel,
    AbstractSoilSource,
    source!

"""
    AbstractSoilBoundaryConditions{FT <: AbstractFloat}

An abstract type for types of boundary conditions, which will include
prescribed functions of space and time as Dirichlet conditions or
Neumann conditions, in addition to other 
convenient soil-specific conditions, like free drainage.
"""
abstract type AbstractSoilBoundaryConditions{FT <: AbstractFloat} end

"""
    AbstractSoilSource{FT <: AbstractFloat}

An abstract type for types of source terms for the soil equations.

In standalone mode, the only supported source type is freezing and 
thawing. ClimaLSM.jl creates additional sources to include as
necessary e.g. root extraction (not available in stand alone mode).
"""
abstract type AbstractSoilSource{FT <: AbstractFloat} end

"""
    AbstractSoilModel{FT} <: AbstractModel{FT} 

The abstract type for all soil models.

Currently, we only have plans to support a RichardsModel, simulating
the flow of liquid water through soil via the Richardson-Richards equation,
 and a fully integrated soil heat
and water model, with phase change.
"""
abstract type AbstractSoilModel{FT} <: ClimaLSM.AbstractModel{FT} end

ClimaLSM.name(::AbstractSoilModel) = :soil

"""
    coordinates(model::AbstractSoilModel)

A extension of the `coordinates` function, which returns the coordinates
of a model domain. 

The coordinates are stored in the model because they are required in
computing the right hand side. 
"""
ClimaLSM.Domains.coordinates(model::AbstractSoilModel) = model.coordinates

"""
   horizontal_components!(dY::ClimaCore.Fields.FieldVector,
                          domain::Column, _...)
Updates dY in place by adding in the tendency terms resulting from
horizontal derivative operators.

In the case of a column domain, there are no horizontal
contributions to the right hand side.
"""
function horizontal_components!(
    dY::ClimaCore.Fields.FieldVector,
    domain::Column,
    _...,
) end

"""
   dss!(dY::ClimaCore.Fields.FieldVector,domain::Column)

Computes the appropriate weighted direct stiffness summation based on
the domain type, updates `dY` in place.

For column domains, no dss is needed.
"""
function dss!(dY::ClimaCore.Fields.FieldVector, domain::Column) end

"""
   dss!(dY::ClimaCore.Fields.FieldVector,domain::HybridBox)

Computes the appropriate weighted direct stiffness summation based on
the domain type, updates `dY` in place.

For the Hybrid box domain, a weighted dss is needed for each variable.
"""
function dss!(dY::ClimaCore.Fields.FieldVector, domain::HybridBox)
    for key in propertynames(dY.soil)
        Spaces.weighted_dss!(getproperty(dY.soil, key))
    end
end

"""
   FluxBC{FT} <: AbstractSoilBoundaryConditions{FT}

A simple concrete type of boundary condition, which enforces
constant fluxes at the top and bottom of the domain.
"""
struct FluxBC{FT} <: AbstractSoilBoundaryConditions{FT}
    top_flux_bc::FT
    bot_flux_bc::FT
end

struct PrecipFreeDrainage{FT} <: AbstractSoilBoundaryConditions{FT}
    top_flux_bc::Function
end


"""
    boundary_fluxes(bc::FluxBC, _...)

A function which returns the correct boundary flux
given the boundary condition type `FluxBC`.

This is a trivial example, but a more complex one would be e.g.
Dirichlet conditions on the state, which then must be converted into
a flux before being applied as a boundary condition.
"""
function boundary_fluxes(bc::FluxBC, _...)
    return bc.top_flux_bc, bc.bot_flux_bc
end
function boundary_fluxes(bc::PrecipFreeDrainage, p, t)
    return bc.top_flux_bc(t), -parent(p.soil.K)[1]
end


"""
     source!(dY::ClimaCore.Fields.FieldVector,
             src::AbstractSoilSource,
             Y::ClimaCore.Fields.FieldVector,
             p::ClimaCore.Fields.FieldVector
             )::ClimaCore.Field.Field

A stub function, which is extended by ClimaLSM.

Once we have the freeze thaw source function, we do not need this
stub anymore.
"""
function source!(
    dY::ClimaCore.Fields.FieldVector,
    src::AbstractSoilSource,
    Y::ClimaCore.Fields.FieldVector,
    p::ClimaCore.Fields.FieldVector,
)::ClimaCore.Field.Field end


include("./soil_hydrology_parameterizations.jl")
include("./soil_heat_parameterizations.jl")
include("./rre.jl")
include("./energy_hydrology.jl")


end

#=
# Domain Tutorial

## Goals of the tutorial
The goal of this is to outline what is currently implemented in ClimaLSM 
and to serve as a software design document 
for future development involving the underlying domains.

## Background

In both the atmosphere and the ocean, all variables are defined at all locations
in the region of interest, or domain.  For example, the air density, temperature, pressure,
and wind speed are defined everywhere in the domain. After choosing a resolution
and discretizing space, the numerical problem is to advance a 
system of ordinary differential
equations, where at each coordinate point a value of 
`ρ`, `T`, `P`, and `u⃗` are solved for at each step. The choice of domain is a question "only"
of geometry: you may be interested in a large eddy simulation (using a box domain), or 
in a global model (where you would need a spherical shell domain
representing the atmosphere or ocean from some depth to a given height).

For land surface
models, each variable is not defined everywhere in space. For example,
the soil water content `θ` is only defined below ground. Snow water equivalent (`S`) is only
defined on the surface itself. Canopy variables are only defined above ground.
 Once we have discretized the land surface region into
a set of points, the numerical problem is to advance a system of ODEs, where
at each coordinate point a different subset of (`θ`, `S`, ...) are solved for.

In other words, different prognostic variables in land surface models exist
in different, overlapping, domains. We need to decide on the geometry of interest (e.g. single column
vs a global simulation), but we also need to specify where each variable of the model is defined.

ClimaLSM Domains were designed with this in mind. Following other Clima models, the domains are defined
so that (1) the user can easily switch geometries, e.g. single column to global model, (2) individual component
models can be run by themselves, using a single domain, and so that (3) different model component domains
 can be combined into a single land model domain, since different land model components can be combined
to run as a single land surface model.

## What is a ClimaLSM Domain?
A domain represents a region of space. In ClimaLSM, domains are simply
structs containing parameters that define these regions - for example an
x-range and y-range that define a plane. In addition, ClimaLSM domains store
some extra information needed to build a ClimaCore function space on the
physical domain, such as the number of elements of a discretization of the
region.

There are two key methods which take a ClimaLSM domain as an argument.

(*) `make_function_space(domain)`: when solving partial differential equations, the spatial 
discretization is tied to a set of basis functions you wish to use to represent the prognostic
variable as a function of space. The nodal points - the locations in space where the variable 
is solved for - are arranged in space in a manner which depends on these basis functions. The method
`make_function_space(domain)` takes in the information stored in the domain (the spatial extent, resolution, etc)
and creates the underlying space. Note that this is only needed for domains of variables satisfying PDEs,
but it is defined for all domains[^1].  


(*) ` coordinates(domain)`: under the hood, this function (1) creates the function space and then (2) uses
the function space to create the coordinate field. This returns the coordinates as a ClimaCore.Fields.Field object.
Depending on the domain, the returned coordinate field will have elements of different names and types. For example,
the SphericalShell domain has coordinates of latitude, longitude, and height, while a Plane domain has coordinates
of x and y, and a Point domain only has a coordinate z_sfc.


## Domain types
All ClimaLSM domains are subtypes of abstract type `ClimaLSM.Domains.AbstractDomain`.
A variety of concrete domain types are supported:

- 0D: `Domains.Point`
- 1D: `Domains.Column`
- 2D: `Domains.Plane`, `Domains.SphericalSurface`
- 3D: `Domains.HybridBox`, `Domains.SphericalShell`.

Single component models (soil, snow, vegetation, canopy airspace...) will be solved on 
a single domain. Which domain is appropriate depends on the model equations (PDE vs ODE) and
on the configuration of interest (single column or 3d). 

Furthermore, these individual model domains are composed when building a multi-component
Land System Model. In this case, the full LSM domain is split into regions:
the subsurface and the surface[^2]. Different land model components are associated
with these different regions - e.g. the soil model is a subsurface model, while vegetation
is a surface model. The above domain types are then paired based on the desired LSM
domain being modeled:

| LSM Domain | Surface Domain | Subsurface Domain |
| :---         |     :---:      |          ---: |
| LSMSingleColumnDomain   | Point    | Column    |
| LSMCartesianBoxDomain    | Plane[^3]       | HybridBox      |
| LSMSphericalShellDomain    | SphericalSurface[^3][^4]       | SphericalShell      |

It is important to note that the horizontal domain used for the surface and subsurface 
domains are identical in LSM simulations. This ensures that we can use the same indexing
of surface and subsurface domains and variables to correctly compute and 
apply boundary fluxes between different component models. Otherwise we would need
to develop additional infrastructure in order to, for example, select the correct subsurface
column corresponding to a particular surface location.



## How variable initialization depends on domains
When a developer creates a model, they need to specify the symbols used for the prognostic variables,
via [`prognostic_vars`](https://clima.github.io/ClimaLSM.jl/dev/APIs/SharedUtilities/#ClimaLSM.prognostic_vars),
and
the types of those variables,
 via [`prognostic_types`](https://clima.github.io/ClimaLSM.jl/dev/APIs/SharedUtilities/#ClimaLSM.prognostic_types).

The [`initialize`](https://clima.github.io/ClimaLSM.jl/dev/APIs/SharedUtilities/#ClimaLSM.initialize)
function (which calls both
[`initialize_prognostic`](https://clima.github.io/ClimaLSM.jl/dev/APIs/SharedUtilities/#ClimaLSM.initialize_prognostic)
 and [`initialize_auxiliary`](https://clima.github.io/ClimaLSM.jl/dev/APIs/SharedUtilities/#ClimaLSM.initialize_auxiliary)
creates the prognostic state vector `Y` (a ClimaCore.Fields.FieldVector). Each field (ClimaCore.Fields.Field) stored within
the field vector corresponds to a prognostic variable (identified with the symbol specified). If the prognostic type for that variable
is a float, the field will be a field of float values (a scalar field)[^5].

How do domains tie into this? The field of prognostic variables corresponds in a 1-1 fashion with the coordinate field created from the
domain. That is, if you coordinate field has points (x,y) in [(1,1), (1,2); (2,1), (2,2)], your prognostic variable field (for variable `θ`,
for example) will be [θ11, θ12; θ21, θ22]. Your variable always has the same spatial resolution as the domain does.


## Future work
Almost all interactions between variables in land surface models are within column - that is, there is only
vertical transport and exchanges. The exception to this is the horizontal flow of water on the surface
and within the soil. The `right hand side` (rhs) functions (the ODE functions) can be split into
"vertical" and "horizontal" pieces. We envision each step of the land surface model simulation to be solved  
in two steps: (1) the vertical rhs evaluations are carried out (and can be parallelized), and (2) the horizontal 
rhs functions are then evaluated (possibly less frequently?) and require communcation between columns.
In this case, right hand side functions will need to be aware of the domain.


[^1]: finite differencing is used in the vertical, and spectral elements are used in the horizontal.

[^2]: a suprasurface region may also be necessary - for example if the canopy airspace model involves PDEs. 

[^3]: We are not sure yet if the surface domain should be associated with a 2d space or if it should be a 3d space but with only one point in the vertical. Both would result in the same number of surface coordinate points, but they would have different underlying representations. The latter may be helpful because we can use the `column` function in the exact same way to extract the same column from the surface and subsurface domains.  In that case, the surface domain would be a HybdridBox (with one vertical element) instead of a Plane, and a SphericalShell (with one vertical element) instead of a Spherical Surface. Regardless, the two domains will share the same instance of the horizontal domain.

[^4] Note that this is not created as a separate object, yet. It is only made in the process of making the SphericalShell domain.

[^5]: We also will support having an array-like type of variable (perhaps an NTuple; this is TBD).  This is to be used for the multi-layer canopy model.
=#

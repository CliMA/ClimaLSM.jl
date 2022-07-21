#=
# Domain Tutorial

## Goals of the tutorial
The goal of this is to outline what is currently implemented in ClimaLSM 
and to serve as a software design document 
for future development involving the underlying domains.

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
but it is defined for all domains[^function_space].  


(*) ` coordinates(domain)`: under the hood, this function (1) creates the function space and then (2) uses
the function space to create the coordinate field. This returns the coordinates as a ClimaCore.Fields.Field object.
Depending on the domain, the returned coordinate field will have elements of different names and types. For example,
the SphericalShell domain has coordinates of latitude, longitude, and height, while a Plane domain has coordinates
of x and y, and a Point domain only has a coordinate z_sfc.


[^function_space]: finite differencing is used in the vertical, and spectral elements are used in the horizontal.

## Examples of domains
All ClimaLSM domains are subtypes of abstract type `ClimaLSM.Domains.AbstractDomain`.
A variety of concrete domain types are supported:

- 0D: `Domains.Point`
- 1D: `Domains.Column`
- 2D: `Domains.Plane`, `Domains.SphericalSurface`
- 3D: `Domains.HybridBox`, `Domains.SphericalShell`.

Single component models (soil, snow, vegetation, canopy airspace...) will be solved on 
a single domain. Which domain is appropriate depends on the model equations (PDE vs ODE) and
on the configuration of interest (single column or 3d). 

Furthermore, these individual model domains can be composed when building a multi-component
Land System Model. In this case, the full LSM domain is split into regions:
the subsurface and the surface[^suprasurfacefootnote]. Different land model components are associated
with these different regions - e.g. the soil model is a subsurface model, while vegetation
is a surface model. The above domain types are then paired based on the desired LSM
domain being modeled:

| LSM Domain | Surface Domain | Subsurface Domain |
| :---         |     :---:      |          ---: |
| LSMSingleColumnDomain   | Point    | Column    |
| LSMCartesianBoxDomain    | Plane[^surfacefootnote]       | HybridBox      |
| LSMSphericalShellDomain    | SphericalSurface[^surfacefootnote][^sphsurffootnote]       | SphericalShell      |

[^suprasurfacefootnote]: a suprasurface region may also be necessary - for example if the
canopy airspace model involves PDEs. 
[^surfacefootnote]: We are not sure yet if the surface domain should be associated with a 2d space or if it should
be a 3d space but with only one point in the vertical. Both would result in the same number of surface
coordinate points, but they would have different underlying representations. The latter may be helpful because
we can use the `column` function in the exact same way to extract the same column from the surface and subsurface domains. 
In that case,
the surface domain would be a HybdridBox (with one vertical element) instead of a Plane, and a SphericalShell
(with one vertical element) instead of a Spherical Surface. Regardless, the two domains will share the same
instance of the horizontal domain.
[^sphsurffootnote] Note that this is not created as a separate object, yet. It is only made in the process of 
making the SphericalShell domain.


## How variable initialization depends on domains
cover initialize, types
Brief reminder of how we initialize - at each coordinate, you get a variable with the type you specify? do we want to discuss this here?
why one then has to specify the domain where interaction variables live
Outline future work on using NTuple as our "array" type

## Future work
Almost all interactions between variables in land surface models are within column - that is, there is only
vertical transport and exchanges. The exception to this is the horizontal flow of water on the surface
and within the soil. The `right hand side` (rhs) functions (the ODE functions) can be split into
"vertical" and "horizontal" pieces. We envision each step of the land surface model simulation to be solved  
in two steps: (1) the vertical rhs evaluations are carried out (and can be parallelized), and (2) the horizontal 
rhs functions are then evaluated (possibly less frequently?) and require communcation between columns.
In this case, right hand side functions will need to be aware of the domain.
=#

using Test
using ClimaCore
if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
using ClimaLSM

using ClimaLSM.Domains: HybridBox, SphericalShell
FT = Float64

@testset "WVector usage in gradient" begin
    # In our RHS, we want to set a boundary condition on flux, which is a gradient
    # of a scalar. This should be a covariant vector. In the case of no topography,
    # the Wvector points in the same direction as the covariant vector.

    # However, they are scaled differently. Check here that supplying a WVector
    # to Operators.DivergenceF2C(top = Operators.SetValue...,) works - i.e.
    # that it is converted internally to a Covariant Vector, and we get the
    # correct output.

    domain = SphericalShell(;
        radius = FT(1.0),
        height = FT(1.0),
        nelements = (1, 2),
        npolynomial = 3,
    )

    coords = ClimaLSM.Domains.coordinates(domain) # center coordinates
    ϕ = (FT(90) .- coords.lat)
    θ = coords.long
    r = coords.z .+ FT(1.0)

    #=
    # For my own sanity, convert lat/lon to usual spherical coords, to cartesian components,
    # and make sure that agrees with CC.
    scalar_x = @.  r*cos(θ*π/180)*sin(ϕ*π/180)
    scalar_y = @.  r*sin(θ*π/180)*sin(ϕ*π/180)
    scalar_z = @. r* cos(ϕ*π/180)
    geom = ClimaCore.Geometry.SphericalGlobalGeometry{FT}(FT(1.0))
    #\vec{r} = r r^
    cartesian = @. ClimaCore.Geometry.CartesianVector(ClimaCore.Geometry.UVWVector(0.0,0.0, r))
    cartesian_points = ClimaCore.Geometry.CartesianPoint.(coords, Ref(geom))
    =#

    scalar = @. r^3 / 1e6
    Δr = 0.5
    upper_bc = ClimaCore.Geometry.Covariant3Vector(2.0 * Δr) # r r^
    lower_bc = ClimaCore.Geometry.Covariant3Vector(1.0 * Δr) # r r^
    gradc2f = ClimaCore.Operators.GradientC2F(
        top = ClimaCore.Operators.SetGradient(upper_bc),
        bottom = ClimaCore.Operators.SetGradient(lower_bc),
    )
    foo = @. gradc2f(scalar)
    divf2c = ClimaCore.Operators.DivergenceF2C()

    div_cov = @. (divf2c(foo))


    upper_bc = ClimaCore.Geometry.Contravariant3Vector(2.0 / Δr) # r r^
    lower_bc = ClimaCore.Geometry.Contravariant3Vector(1.0 / Δr) # r r^
    gradc2f = ClimaCore.Operators.GradientC2F(
        top = ClimaCore.Operators.SetGradient(upper_bc),
        bottom = ClimaCore.Operators.SetGradient(lower_bc),
    )
    foo = @. gradc2f(scalar)
    divf2c = ClimaCore.Operators.DivergenceF2C()

    div_con = @. (divf2c(foo))

    upper_bc = ClimaCore.Geometry.WVector(2.0) # r r^
    lower_bc = ClimaCore.Geometry.WVector(1.0) # r r^
    gradc2f = ClimaCore.Operators.GradientC2F(
        top = ClimaCore.Operators.SetGradient(upper_bc),
        bottom = ClimaCore.Operators.SetGradient(lower_bc),
    )
    foo = @. gradc2f(scalar)
    divf2c = ClimaCore.Operators.DivergenceF2C()

    div_w = @. (divf2c(foo))
    @test parent(div_w) == parent(div_cov)
    @test parent(div_cov) == parent(div_con)


    # Unclear why this does not work!
    # upper_bc =ClimaCore.Geometry.Covariant3Vector(2.0*Δr) # r r^
    # lower_bc =ClimaCore.Geometry.Covariant3Vector(1.0*Δr) # r r^
    # gradc2f = ClimaCore.Operators.GradientC2F()
    # divf2c = ClimaCore.Operators.DivergenceF2C(
    #     top = ClimaCore.Operators.SetValue(upper_bc),
    #     bottom = ClimaCore.Operators.SetValue(lower_bc),
    # )

    #div_rpt = unique(parent(@. (divf2c(gradc2f(scalar)))))
end

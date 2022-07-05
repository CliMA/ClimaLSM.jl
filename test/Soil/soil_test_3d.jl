using Test
using UnPack
using Statistics
using ClimaCore
import CLIMAParameters as CP

if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
using ClimaLSM
using ClimaLSM.Domains: HybridBox, SphericalShell
using ClimaLSM.Soil
import ClimaLSM
include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))

FT = Float64

@testset "Soil horizontal operators unit tests" begin
    ν = FT(0.495)
    K_sat = FT(1)#0.0443 / 3600 / 100); # m/s
    S_s = FT(1)#e-3); #inverse meters
    vg_n = FT(2.0)
    vg_α = FT(2.6) # inverse meters
    vg_m = FT(1) - FT(1) / vg_n
    θ_r = FT(0)
    zmax = FT(0)
    zmin = FT(-1)
    xmax = FT(1.0)
    soil_domain = HybridBox(;
        xlim = (0.0, xmax),
        ylim = (0.0, 1.0),
        zlim = (zmin, zmax),
        nelements = (500, 2, 500),
        npolynomial = 3,
    )
    top_flux_bc = FT(0.0)
    bot_flux_bc = FT(0.0)
    sources = ()
    boundary_fluxes = FluxBC{FT}(top_flux_bc, bot_flux_bc)
    params = Soil.RichardsParameters{FT}(ν, vg_α, vg_n, vg_m, K_sat, S_s, θ_r)

    soil = Soil.RichardsModel{FT}(;
        parameters = params,
        domain = soil_domain,
        boundary_conditions = boundary_fluxes,
        sources = sources,
    )

    # sinusoidal water table
    function init_soil!(Ysoil, x, z, params)
        function hydrostatic_profile(
            x::FT,
            z::FT,
            params::RichardsParameters{FT},
        ) where {FT}
            @unpack ν, vg_α, vg_n, vg_m, θ_r, S_s = params
            z_∇ = FT(zmin / 2.0 + (zmax - zmin) / 10.0 * sin(π * 2 * x / xmax))
            if z > z_∇
                S = FT((FT(1) + (vg_α * (z - z_∇))^vg_n)^(-vg_m))
                ϑ_l = S * (ν - θ_r) + θ_r
            else
                ϑ_l = -S_s * (z - z_∇) + ν
            end
            return FT(ϑ_l)
        end
        Ysoil.soil.ϑ_l .= hydrostatic_profile.(x, z, Ref(params))
    end
    soil_ode! = make_ode_function(soil)
    Y, p, coords = initialize(soil)
    init_soil!(Y, coords.x, coords.z, soil.parameters)
    dY = similar(Y)
    soil_ode!(dY, Y, p, 0.0)

    function dθdx(x, z)
        z_∇ = zmin / 2.0 + (zmax - zmin) / 10.0 * sin(π * 2 * x / xmax)
        dz∇dx = (zmax - zmin) / 10.0 * cos(π * 2 * x / xmax) * (2 * π / xmax)
        if z > z_∇
            dSdz_∇ =
                -(vg_m) *
                (1 + (vg_α * (z - z_∇))^vg_n)^(-vg_m - 1) *
                vg_n *
                (vg_α * (z - z_∇))^(vg_n - 1) *
                (-vg_α)
            dSdx = dSdz_∇ * dz∇dx
            return (ν - θ_r) * dSdx
        else
            return S_s * dz∇dx
        end
    end

    function d2θdx2(x, z)
        z_∇ = zmin / 2.0 + (zmax - zmin) / 10.0 * sin(π * 2 * x / xmax)
        dz∇dx = (zmax - zmin) / 10.0 * cos(π * 2 * x / xmax) * (2π / xmax)
        d2z∇dx2 =
            -(zmax - zmin) / 10.0 * sin(π * 2 * x / xmax) * (2π / xmax)^2.0
        if z > z_∇
            dSdz_∇ =
                -(vg_m) *
                (1 + (vg_α * (z - z_∇))^vg_n)^(-vg_m - 1) *
                vg_n *
                (vg_α * (z - z_∇))^(vg_n - 1) *
                (-vg_α)

            d2Sdz∇2 =
                -vg_m *
                vg_n *
                vg_α^2.0 *
                (
                    (vg_n - 1) *
                    (vg_α * (z - z_∇))^(vg_n - 2) *
                    (1 + (vg_α * (z - z_∇))^vg_n)^(-vg_m - 1) -
                    (vg_m + 1) *
                    vg_n *
                    (vg_α * (z - z_∇))^(2.0 * (vg_n - 1)) *
                    (1 + (vg_α * (z - z_∇))^vg_n)^(-vg_m - 2)
                )
            d2Sdx2 = dz∇dx^2.0 * d2Sdz∇2 + dSdz_∇ * d2z∇dx2
            return (ν - θ_r) * d2Sdx2
        else
            return S_s * d2z∇dx2
        end
    end

    function dKdθ(θ)
        S = (θ - θ_r) / (ν - θ_r)
        if S < 1
            f = 1.0 - (1.0 - S^(1.0 / vg_m))^vg_m
            f1 = f^2.0 / 2.0 / S^0.5
            f2 = 2 * S^(1 / vg_m - 1 / 2) * f / (1 - S^(1 / vg_m))^(1.0 - vg_m)
            return (f1 + f2) * K_sat / (ν - θ_r)
        else
            return 0.0
        end

    end

    function dψdθ(θ)
        S = (θ - θ_r) / (ν - θ_r)
        if S < 1.0
            return 1.0 / (vg_α * vg_m * vg_n) / (ν - θ_r) *
                   (S^(-1 / vg_m) - 1)^(1 / vg_n - 1) *
                   S^(-1 / vg_m - 1)
        else
            return 1.0 / S_s
        end
    end


    function d2ψdθ2(θ)
        S = (θ - θ_r) / (ν - θ_r)
        if S < 1.0
            return 1.0 / (vg_α * vg_m * vg_n) / (ν - θ_r)^2.0 * (
                S^(-2.0 / vg_m - 2.0) *
                (-1 / vg_m) *
                (1 / vg_n - 1) *
                (S^(-1 / vg_m) - 1)^(1 / vg_n - 2) +
                (S^(-1 / vg_m) - 1)^(1 / vg_n - 1) *
                (-1 / vg_m - 1) *
                S^(-1 / vg_m - 2)
            )
        else
            return 0.0
        end
    end

    function K(θ)
        S = (θ - θ_r) / (ν - θ_r)
        if S < 1
            f = 1.0 - (1.0 - S^(1.0 / vg_m))^vg_m
            return K_sat * f^2.0 * sqrt(S)
        else
            return K_sat
        end

    end


    X = coords.x
    Z = coords.z
    θ = Y.soil.ϑ_l

    #5.517176201418359e-9 max error with horizontal terms off

    # Test in unsaturated zone
    for N in [150, 350]
        myXslice = parent(X)[N, 1, 1, 1, :]
        myθslice = parent(Y.soil.ϑ_l)[N, 1, 1, 1, :]
        mydYslice = parent(dY.soil.ϑ_l)[N, 1, 1, 1, :]
        local_p = sortperm(myXslice)

        Xsort = myXslice[local_p]
        θsort = myθslice[local_p]
        dYsort = mydYslice[local_p]
        unique_indices = unique(i -> Xsort[i], 1:length(Xsort))
        XX = Xsort[unique_indices]
        θX = θsort[unique_indices]
        myZ = unique(parent(Z)[N, 1, 1, 1, :])[1]
        dYX = dYsort[unique_indices]
        expected = @. (
            dKdθ(θX) * dψdθ(θX) * dθdx(XX, myZ)^2.0 +
            K(θX) * dψdθ(θX) * d2θdx2(XX, myZ) +
            K(θX) * dθdx(XX, myZ)^2.0 * d2ψdθ2(θX)
        )
        @test maximum(abs.(dYX .- expected)) < 1e-8
    end
end


@testset "Soil energy+hydrology horizontal operators" begin
    earth_param_set = create_lsm_parameters(FT)
    ν = FT(0.44)
    K_sat = FT(29.7 / 3600 / 100) # m/s
    S_s = FT(1e-3) #inverse meters
    vg_n = FT(2.68)
    vg_α = FT(14.5) # inverse meters
    vg_m = FT(1) - FT(1) / vg_n
    θ_r = FT(0.045)
    ν_ss_om = FT(0.0)
    ν_ss_quartz = FT(1.0)
    ν_ss_gravel = FT(0.0)
    κ_minerals = FT(2.5)
    κ_om = FT(0.25)
    κ_quartz = FT(8.0)
    κ_air = FT(0.025)
    κ_ice = FT(2.21)
    κ_liq = FT(0.57)
    ρp = FT(2.66 / 1e3 * 1e6)
    ρc_ds = FT(2e6)
    κ_solid = Soil.κ_solid(ν_ss_om, ν_ss_quartz, κ_om, κ_quartz, κ_minerals)
    κ_dry = Soil.κ_dry(ρp, ν, κ_solid, κ_air)
    κ_sat_frozen = Soil.κ_sat_frozen(κ_solid, ν, κ_ice)
    κ_sat_unfrozen = Soil.κ_sat_unfrozen(κ_solid, ν, κ_liq)
    parameters = Soil.EnergyHydrologyParameters(
        κ_dry = κ_dry,
        κ_sat_frozen = κ_sat_frozen,
        κ_sat_unfrozen = κ_sat_unfrozen,
        ρc_ds = ρc_ds,
        ν = ν,
        ν_ss_om = ν_ss_om,
        ν_ss_quartz = ν_ss_quartz,
        ν_ss_gravel = ν_ss_gravel,
        vg_α = vg_α,
        vg_n = vg_n,
        K_sat = K_sat,
        S_s = S_s,
        θ_r = θ_r,
        earth_param_set = earth_param_set,
    )
    zmax = FT(0)
    zmin = FT(-1)
    xmax = FT(1.0)
    soil_domain = HybridBox(;
        xlim = (0.0, xmax),
        ylim = (0.0, 1.0),
        zlim = (zmin, zmax),
        nelements = (500, 2, 500),
        npolynomial = 3,
    )
    top_flux_bc = FT(0.0)
    bot_flux_bc = FT(0.0)

    sources = ()
    boundary_fluxes = Soil.FluxBC{FT}(top_flux_bc, bot_flux_bc)
    soil = Soil.EnergyHydrology{FT}(;
        parameters = parameters,
        domain = soil_domain,
        rre_boundary_conditions = boundary_fluxes,
        heat_boundary_conditions = boundary_fluxes,
        sources = sources,
    )

    Y, p, coords = initialize(soil)

    # specify ICs
    function init_soil!(Ysoil, z, params)
        function hydrostatic_profile(
            z::FT,
            params::EnergyHydrologyParameters,
        ) where {FT}
            @unpack ν, vg_α, vg_n, vg_m, θ_r, S_s = params
            z_∇ = FT(zmin / 2.0)
            if z > z_∇
                S = FT((FT(1) + (vg_α * (z - z_∇))^vg_n)^(-vg_m))
                ϑ_l = S * (ν - θ_r) + θ_r
            else
                ϑ_l = -S_s * (z - z_∇) + ν
            end
            return FT(ϑ_l)
        end
        Ysoil.soil.ϑ_l .= hydrostatic_profile.(z, Ref(params))
        Ysoil.soil.θ_i .= ClimaCore.Fields.zeros(FT, axes(Ysoil.soil.θ_i))
    end
    init_soil!(Y, coords.z, soil.parameters)


    θ_l = Soil.volumetric_liquid_fraction.(Y.soil.ϑ_l, ν)
    volumetric_heat_capacity =
        Soil.volumetric_heat_capacity.(θ_l, Y.soil.θ_i, Ref(parameters))
    Y.soil.ρe_int .=
        ClimaCore.Fields.zeros(FT, axes(Y.soil.ρe_int)) .+
        volumetric_internal_energy.(
            0.0,
            volumetric_heat_capacity,
            280.0,
            Ref(parameters),
        )



    soil_ode! = make_ode_function(soil)
    dY = similar(Y)
    soil_ode!(dY, Y, p, 0.0)
    ## dY should be zero. Look at dY/Y. 
    @test maximum(abs.(parent(dY.soil.ϑ_l))) / ν < 5e-12
    @test maximum(abs.(parent(dY.soil.θ_i))) / ν < 5e-12
    @test maximum(abs.(parent(dY.soil.ρe_int))) ./ 2.052e7 < 5e-12
end


@testset "Soil hydrology RHS on sphere" begin
    ν = FT(0.44)
    K_sat = FT(1.0)
    S_s = FT(1e-3) #inverse meters
    vg_n = FT(2.68)
    vg_α = FT(14.5) # inverse meters
    vg_m = FT(1) - FT(1) / vg_n
    θ_r = FT(0.045)

    parameters = Soil.RichardsParameters(
        ν = ν,
        vg_α = vg_α,
        vg_n = vg_n,
        vg_m = vg_m,
        K_sat = K_sat,
        S_s = S_s,
        θ_r = θ_r,
    )

    soil_domain = SphericalShell(;
        radius = FT(100.0),
        height = FT(1.0),
        nelements = (1, 30),
        npolynomial = 3,
    )
    top_flux_bc = FT(K_sat)
    bot_flux_bc = FT(K_sat)

    sources = ()
    boundary_fluxes = Soil.FluxBC{FT}(top_flux_bc, bot_flux_bc)
    soil = Soil.RichardsModel{FT}(;
        parameters = parameters,
        domain = soil_domain,
        boundary_conditions = boundary_fluxes,
        sources = sources,
    )

    Y, p, coords = initialize(soil)

    # specify ICs
    function init_soil!(Ysoil, z, params)
        function hydrostatic_profile(
            z::FT,
            params::RichardsParameters,
        ) where {FT}
            @unpack ν, vg_α, vg_n, vg_m, θ_r, S_s = params
            z_∇ = FT(0.5)
            if z > z_∇
                S = FT((FT(1) + (vg_α * (z - z_∇))^vg_n)^(-vg_m))
                ϑ_l = S * (ν - θ_r) + θ_r
            else
                ϑ_l = -S_s * (z - z_∇) + ν
            end
            return FT(ϑ_l)
        end
        Ysoil.soil.ϑ_l .= hydrostatic_profile.(z, Ref(params))
    end
    init_soil!(Y, coords.z, soil.parameters)
    soil_ode! = make_ode_function(soil)
    dY = similar(Y)
    soil_ode!(dY, Y, p, 0.0)
    ## dY should be zero except at the boundary, where it should be ±30.
    @test (mean(unique(parent(ClimaCore.level(dY.soil.ϑ_l, 1)))) - 30.0) / ν <
          1e-11
    @test (mean(unique(parent(ClimaCore.level(dY.soil.ϑ_l, 30)))) + 30.0) / ν <
          1e-11
    @test mean([
        maximum(abs.(unique(parent(ClimaCore.level(dY.soil.ϑ_l, k))))) for
        k in 2:29
    ]) / ν < 1e-11

end

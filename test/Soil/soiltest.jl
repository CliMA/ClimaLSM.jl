using Test
using UnPack
using Statistics
using DiffEqCallbacks
using OrdinaryDiffEq: ODEProblem, solve, Euler
using ClimaCore
import CLIMAParameters as CP

if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
using ClimaLSM
using ClimaLSM.Domains: Column
using ClimaLSM.Soil

import ClimaLSM
include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))

FT = Float64


@testset "Richards equation" begin

    saved_values = SavedValues(FT, ClimaCore.Fields.FieldVector)
    ν = FT(0.495)
    K_sat = FT(0.0443 / 3600 / 100) # m/s
    S_s = FT(1e-3) #inverse meters
    vg_n = FT(2.0)
    vg_α = FT(2.6) # inverse meters
    vg_m = FT(1) - FT(1) / vg_n
    θ_r = FT(0)
    zmax = FT(0)
    zmin = FT(-10)
    nelems = 50

    soil_domain = Column(; zlim = (zmin, zmax), nelements = nelems)
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

    Y, p, coords = initialize(soil)

    # specify ICs
    function init_soil!(Ysoil, z, params)
        function hydrostatic_profile(
            z::FT,
            params::RichardsParameters{FT},
        ) where {FT}
            @unpack ν, vg_α, vg_n, vg_m, θ_r = params
            #unsaturated zone only, assumes water table starts at z_∇
            z_∇ = FT(-10)# matches zmin
            S = FT((FT(1) + (vg_α * (z - z_∇))^vg_n)^(-vg_m))
            ϑ_l = S * (ν - θ_r) + θ_r
            return FT(ϑ_l)
        end
        Ysoil.soil.ϑ_l .= hydrostatic_profile.(z, Ref(params))
    end

    init_soil!(Y, coords.z, soil.parameters)

    soil_ode! = make_ode_function(soil)

    t0 = FT(0)
    tf = FT(60)
    dt = FT(1)
    cb = SavingCallback((u, t, integrator) -> copy(integrator.p), saved_values)
    prob = ODEProblem(soil_ode!, Y, (t0, tf), p)
    sol = solve(prob, Euler(); dt = dt, callback = cb)

    @test sum(parent(sol.u[end]) .== parent(Y.soil.ϑ_l)) == nelems
    # should be hydrostatic equilibrium at every layer, at each step:
    @test mean(
        sum([
            parent(saved_values.saveval[k].soil.ψ .+ coords.z)[:] .+ 10.0 for
            k in 2:1:50
        ]),
    ) < 1e-10

end


@testset "Soil Energy and Water RHS unit tests" begin
    earth_param_set = create_lsm_parameters(FT)
    ν = FT(0.495)
    K_sat = FT(0.0443 / 3600 / 100) # m/s
    S_s = FT(1e-3) #inverse meters
    vg_n = FT(2.0)
    vg_α = FT(2.6) # inverse meters
    vg_m = FT(1) - FT(1) / vg_n
    θ_r = FT(0.1)
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
    ρc_ds = FT(2e6 * (1.0 - ν))
    κ_solid = Soil.κ_solid(ν_ss_om, ν_ss_quartz, κ_om, κ_quartz, κ_minerals)
    κ_dry = Soil.κ_dry(ρp, ν, κ_solid, κ_air)
    κ_sat_frozen = Soil.κ_sat_frozen(κ_solid, ν, κ_ice)
    κ_sat_unfrozen = Soil.κ_sat_unfrozen(κ_solid, ν, κ_liq)
    zmax = FT(0)
    zmin = FT(-1)
    nelems = 200
    Δz = 0.005
    soil_domain = Column(; zlim = (zmin, zmax), nelements = nelems)
    top_flux_bc = FT(0.0)
    bot_flux_bc = FT(0.0)

    sources = ()
    boundary_fluxes = Soil.FluxBC{FT}(top_flux_bc, bot_flux_bc)


    ###
    hyd_off_en_on = Soil.EnergyHydrologyParameters(;
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
        K_sat = 0.0,
        S_s = S_s,
        θ_r = θ_r,
        earth_param_set = earth_param_set,
    )
    soil_heat_on = Soil.EnergyHydrology{FT}(;
        parameters = hyd_off_en_on,
        domain = soil_domain,
        rre_boundary_conditions = boundary_fluxes,
        heat_boundary_conditions = boundary_fluxes,
        sources = sources,
    )

    Y, p, coords = initialize(soil_heat_on)

    # specify ICs
    function init_soil_heat!(Ysoil, z, params)
        ν = params.ν
        Ysoil.soil.ϑ_l .= ν / 2.0
        Ysoil.soil.θ_i .= 0.0
        T = 280.0 .+ 0.5 .* (z .+ 0.5) .^ 2.0
        ρc_s = Soil.volumetric_heat_capacity(ν / 2.0, 0.0, params)
        Ysoil.soil.ρe_int =
            Soil.volumetric_internal_energy.(0.0, ρc_s, T, Ref(params))
    end

    init_soil_heat!(Y, coords.z, soil_heat_on.parameters)
    soil_ode! = make_ode_function(soil_heat_on)
    dY = similar(Y)
    soil_ode!(dY, Y, p, 0.0)
    F_face = 0.0
    κ = parent(p.soil.κ)
    F_below = -0.5 * (κ[end] + κ[end - 1]) * (-Δz + 0.5)
    dY_top = -(F_face - F_below) / Δz
    F_top = -0.5 * (κ[2] + κ[1]) * ((-1.0 + Δz) + 0.5)
    dY_bot = -(F_top - F_face) / Δz
    expected = parent(p.soil.κ)
    expected[1] = dY_bot
    expected[end] = dY_top
    @test mean(abs.(expected .- parent(dY.soil.ρe_int))) /
          median(parent(Y.soil.ρe_int)) < 1e-13
    @test maximum(abs.(parent(dY.soil.ϑ_l))) == 0.0
    @test maximum(abs.(parent(dY.soil.θ_i))) == 0.0

    ###
    hyd_on_en_off = Soil.EnergyHydrologyParameters(
        κ_dry = 0.0,
        κ_sat_frozen = 0.0,
        κ_sat_unfrozen = 0.0,
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
    soil_water_on = Soil.EnergyHydrology{FT}(;
        parameters = hyd_on_en_off,
        domain = soil_domain,
        rre_boundary_conditions = boundary_fluxes,
        heat_boundary_conditions = boundary_fluxes,
        sources = sources,
    )

    Y, p, coords = initialize(soil_water_on)

    # specify ICs
    function init_soil_water!(Ysoil, z, params)
        ν = params.ν
        Ysoil.soil.ϑ_l .= ν / 2.0 .+ ν / 4.0 .* (z .+ 0.5) .^ 2.0
        Ysoil.soil.θ_i .= 0.0
        ρc_s = Soil.volumetric_heat_capacity.(Y.soil.ϑ_l, 0.0, Ref(params))

        Ysoil.soil.ρe_int =
            volumetric_internal_energy.(0.0, ρc_s, 288.0, Ref(params))
    end

    init_soil_water!(Y, coords.z, soil_water_on.parameters)
    soil_ode! = make_ode_function(soil_water_on)
    dY = similar(Y)
    soil_ode!(dY, Y, p, 0.0)
    function dKdθ(θ::FT)::FT where {FT}
        S = (θ - θ_r) / (ν - θ_r)
        if S < 1
            f::FT = 1.0 - (1.0 - S^(1.0 / vg_m))^vg_m
            f1::FT = f^2.0 / 2.0 / S^0.5
            f2::FT =
                2 * S^(1 / vg_m - 1 / 2) * f / (1 - S^(1 / vg_m))^(1.0 - vg_m)
            return (f1 + f2) * K_sat / (ν - θ_r)
        else
            return 0.0
        end

    end

    function dψdθ(θ::FT)::FT where {FT}
        S = (θ - θ_r) / (ν - θ_r)
        if S < 1.0
            return 1.0 / (vg_α * vg_m * vg_n) / (ν - θ_r) *
                   (S^(-1 / vg_m) - 1)^(1 / vg_n - 1) *
                   S^(-1 / vg_m - 1)
        else
            return 1.0 / S_s
        end
    end


    function d2ψdθ2(θ::FT)::FT where {FT}
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

    function K(θ::FT)::FT where {FT}
        S = (θ - θ_r) / (ν - θ_r)
        if S < 1
            f = 1.0 - (1.0 - S^(1.0 / vg_m))^vg_m
            return K_sat * f^2.0 * sqrt(S)
        else
            return K_sat
        end

    end

    function dθdz(z::FT)::FT where {FT}
        return ν / 2.0 * (z + 0.5)
    end


    function d2θdz2(z::FT)::FT where {FT}
        return ν / 2.0
    end

    θ = parent(Y.soil.ϑ_l)# on the center
    θ_face = 0.5 * (θ[2:end] + θ[1:(end - 1)])
    Z = parent(coords.z)
    Z_face = 0.5 * (Z[2:end] + Z[1:(end - 1)])
    K_face = 0.5 .* (K.(θ[2:end]) .+ K.(θ[1:(end - 1)]))
    flux_interior = (@. -K_face * (1.0 + dψdθ(θ_face) * dθdz(Z_face)))
    flux = vcat([0.0], flux_interior, [0.0])

    expected = -(flux[2:end] - flux[1:(end - 1)]) ./ Δz
    @test mean(abs.(expected .- parent(dY.soil.ϑ_l))) / ν < 1e-13
    @test maximum(abs.(parent(dY.soil.θ_i))) == 0.0

    ρe_int_l = parent(
        Soil.volumetric_internal_energy_liq.(
            p.soil.T,
            Ref(soil_water_on.parameters),
        ),
    )
    ρe_int_l_face = 0.5 * (ρe_int_l[2:end] + ρe_int_l[1:(end - 1)])
    flux_interior =
        (@. -K_face * ρe_int_l_face * (1.0 + dψdθ(θ_face) * dθdz(Z_face)))
    flux = vcat([0.0], flux_interior, [0.0])
    expected = -(flux[2:end] - flux[1:(end - 1)]) ./ Δz

    @test mean(abs.(expected .- parent(dY.soil.ρe_int))) /
          median(parent(Y.soil.ρe_int)) < 1e-13

    ###

    hyd_off_en_off = Soil.EnergyHydrologyParameters(;
        κ_dry = 0.0,
        κ_sat_frozen = 0.0,
        κ_sat_unfrozen = 0.0,
        ρc_ds = ρc_ds,
        ν = ν,
        ν_ss_om = ν_ss_om,
        ν_ss_quartz = ν_ss_quartz,
        ν_ss_gravel = ν_ss_gravel,
        vg_α = vg_α,
        vg_n = vg_n,
        K_sat = 0.0,
        S_s = S_s,
        θ_r = θ_r,
        earth_param_set = earth_param_set,
    )

    soil_both_off = Soil.EnergyHydrology{FT}(;
        parameters = hyd_off_en_off,
        domain = soil_domain,
        rre_boundary_conditions = boundary_fluxes,
        heat_boundary_conditions = boundary_fluxes,
        sources = sources,
    )

    Y, p, coords = initialize(soil_both_off)

    # specify ICs
    function init_soil_off!(Ysoil, z, params)
        ν = params.ν
        Ysoil.soil.ϑ_l .= ν / 2.0
        Ysoil.soil.θ_i .= 0.0
        T = 280.0 .+ 0.5 .* (z .+ 0.5) .^ 2.0
        ρc_s = Soil.volumetric_heat_capacity(ν / 2.0, 0.0, params)
        Ysoil.soil.ρe_int =
            Soil.volumetric_internal_energy.(0.0, ρc_s, T, Ref(params))
    end

    init_soil_off!(Y, coords.z, soil_both_off.parameters)
    soil_ode! = make_ode_function(soil_both_off)
    dY = similar(Y)
    soil_ode!(dY, Y, p, 0.0)
    @test maximum(abs.(parent(dY.soil.ρe_int))) == 0.0
    @test maximum(abs.(parent(dY.soil.ϑ_l))) == 0.0
    @test maximum(abs.(parent(dY.soil.θ_i))) == 0.0


    ### Test with both on

    hyd_on_en_on = Soil.EnergyHydrologyParameters(;
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

    soil_both_on = Soil.EnergyHydrology{FT}(;
        parameters = hyd_on_en_on,
        domain = soil_domain,
        rre_boundary_conditions = boundary_fluxes,
        heat_boundary_conditions = boundary_fluxes,
        sources = sources,
    )

    Y, p, coords = initialize(soil_both_on)

    # specify ICs
    function init_soil_on!(Ysoil, z, params)
        ν = params.ν
        Ysoil.soil.ϑ_l .= ν / 2.0 .+ ν / 4.0 .* (z .+ 0.5) .^ 2.0
        Ysoil.soil.θ_i .= 0.0
        T = 280.0 .+ 0.5 .* (z .+ 0.5) .^ 2.0
        ρc_s = Soil.volumetric_heat_capacity.(Y.soil.ϑ_l, 0.0, Ref(params))
        Ysoil.soil.ρe_int =
            volumetric_internal_energy.(0.0, ρc_s, T, Ref(params))
    end

    init_soil_on!(Y, coords.z, soil_both_on.parameters)
    soil_ode! = make_ode_function(soil_both_on)
    dY = similar(Y)
    soil_ode!(dY, Y, p, 0.0)
    @test maximum(abs.(parent(dY.soil.θ_i))) == 0.0

    θ = parent(Y.soil.ϑ_l)# on the center
    θ_face = 0.5 * (θ[2:end] + θ[1:(end - 1)])
    Z = parent(coords.z)
    Z_face = 0.5 * (Z[2:end] + Z[1:(end - 1)])
    K_face = 0.5 .* (K.(θ[2:end]) .+ K.(θ[1:(end - 1)]))
    vf = parent(
        Soil.viscosity_factor.(p.soil.T, hyd_on_en_on.γ, hyd_on_en_on.γT_ref),
    )
    vf_face = 0.5 .* (vf[2:end] .+ vf[1:(end - 1)])

    flux_interior = (@. -K_face * vf_face * (1.0 + dψdθ(θ_face) * dθdz(Z_face)))
    flux = vcat([0.0], flux_interior, [0.0])

    expected = -(flux[2:end] - flux[1:(end - 1)]) ./ Δz
    @test mean(abs.(expected .- parent(dY.soil.ϑ_l))) / ν < 1e-13

    ρe_int_l = parent(
        Soil.volumetric_internal_energy_liq.(
            p.soil.T,
            Ref(soil_both_on.parameters),
        ),
    )
    ρe_int_l_face = 0.5 * (ρe_int_l[2:end] + ρe_int_l[1:(end - 1)])
    κ = parent(p.soil.κ)
    κ_face = 0.5 * (κ[2:end] + κ[1:(end - 1)])


    flux_interior = (@. -K_face *
        vf_face *
        ρe_int_l_face *
        (1.0 + dψdθ(θ_face) * dθdz(Z_face)))
    flux_one = vcat([0.0], flux_interior, [0.0])
    flux_interior = (@. -κ_face * (Z_face + 0.5))
    flux_two = vcat([0.0], flux_interior, [0.0])
    flux = flux_one .+ flux_two
    expected = -(flux[2:end] - flux[1:(end - 1)]) ./ Δz

    @test mean(abs.(expected .- parent(dY.soil.ρe_int))) /
          median(parent(Y.soil.ρe_int)) < 1e-13

end

using Test
using DiffEqCallbacks
using UnPack
using OrdinaryDiffEq: ODEProblem, solve, Euler
using ClimaCore
import CLIMAParameters as CP

if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end

using ClimaLSM
using ClimaLSM.Domains: Column, RootDomain
using ClimaLSM.Soil
using ClimaLSM.Roots
import ClimaLSM
include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))

FT = Float64
@testset "Root soil LSM integration test" begin
    saved_values = SavedValues(FT, ClimaCore.Fields.FieldVector)
    earth_param_set = create_lsm_parameters(FT)
    a_root = FT(13192)
    a_stem = FT(515.5605)
    b_root = FT(2.1079)
    b_stem = FT(0.9631)
    size_reservoir_leaf_moles = FT(16766.2790)
    size_reservoir_stem_moles = FT(11000.8837)
    K_max_root_moles = FT(12.9216)
    K_max_stem_moles = FT(3.4415)
    z_leaf = FT(12) # height of leaf
    # currently hardcoded to match the soil coordinates. this has to
    # be fixed eventually.
    z_root_depths = -Array(1:1:20.0) ./ 20.0 * 3.0 .+ 0.15 / 2.0
    z_bottom_stem = FT(0.0)
    roots_domain = RootDomain{FT}(z_root_depths, [z_bottom_stem, z_leaf])
    roots_ps = Roots.RootsParameters{FT, typeof(earth_param_set)}(
        a_root,
        b_root,
        a_stem,
        b_stem,
        size_reservoir_stem_moles,
        size_reservoir_leaf_moles,
        K_max_root_moles,
        K_max_stem_moles,
        earth_param_set,
    )

    zmin = FT(-3.0)
    zmax = FT(0.0)
    nelements = 20
    soil_domain = Column(; zlim = (zmin, zmax), nelements = nelements)
    ν = FT(0.495)
    K_sat = FT(0.0443 / 3600 / 100) # m/s
    S_s = FT(1e-3) #inverse meters
    vg_n = FT(2.0)
    vg_α = FT(2.6) # inverse meters
    vg_m = FT(1) - FT(1) / vg_n
    θ_r = FT(0)
    soil_ps = Soil.RichardsParameters{FT}(ν, vg_α, vg_n, vg_m, K_sat, S_s, θ_r)
    soil_args = (domain = soil_domain, parameters = soil_ps)
    root_args = (domain = roots_domain, parameters = roots_ps)
    land = RootSoilModel{FT}(;
        soil_model_type = Soil.RichardsModel{FT},
        soil_args = soil_args,
        vegetation_model_type = Roots.RootsModel{FT},
        vegetation_args = root_args,
    )
    Y, p, coords = initialize(land)
    # specify ICs
    function init_soil!(Ysoil, z, params)
        function hydrostatic_profile(
            z::FT,
            params::RichardsParameters{FT},
        ) where {FT}
            @unpack ν, vg_α, vg_n, vg_m, θ_r = params
            #unsaturated zone only, assumes water table starts at z_∇
            z_∇ = FT(-3)# matches zmin
            S = FT((FT(1) + (vg_α * (z - z_∇))^vg_n)^(-vg_m))
            ϑ_l = S * (ν - θ_r) + θ_r
            return FT(ϑ_l)
        end
        Ysoil.soil.ϑ_l .= hydrostatic_profile.(z, Ref(params))
    end
    init_soil!(Y, coords.subsurface.z, land.soil.parameters)

    ## soil is at total ψ+z = -3.0 #m
    ## Want ρgΨ_plant = ρg(-3) - ρg z_plant & convert to MPa
    # we should standardize the units! and not ahve to convert every time.
    # convert parameters once up front and then not each RHS
    p_stem_ini = (-3.0 - z_bottom_stem) * 9.8 * 1000.0 / 1000000.0
    p_leaf_ini = (-3.0 - z_leaf) * 9.8 * 1000.0 / 1000000.0

    theta_stem_0 = p_to_theta(p_stem_ini)
    theta_leaf_0 = p_to_theta(p_leaf_ini)
    y1_0 = FT(theta_stem_0 * size_reservoir_stem_moles)
    y2_0 = FT(theta_leaf_0 * size_reservoir_leaf_moles)
    y0 = [y1_0, y2_0]
    Y.vegetation.rwc .= y0

    ode! = make_ode_function(land)
    t0 = FT(0)
    tf = FT(2)
    dt = FT(1)
    cb = SavingCallback((u, t, integrator) -> integrator.p, saved_values)
    prob = ODEProblem(ode!, Y, (t0, tf), p)
    sol = solve(prob, Euler(), dt = dt, callback = cb)
    #Currently just testing to make sure it runs, but need to have a better test suite.
end

using Test
using DifferentialEquations
using UnPack
using OrdinaryDiffEq: ODEProblem, solve, RK4
using ClimaCore

if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end

using ClimaLSM
using ClimaLSM.Domains: Column, RootDomain
using ClimaLSM.Soil
using ClimaLSM.Roots

FT = Float64
@testset "Root soil LSM integration test" begin
    saved_values = SavedValues(FT, ClimaCore.Fields.FieldVector)
    a_root = FT(0.1)
    a_stem = a_root
    b_root = FT(0.17 / 1e6) # Inverse Pa
    b_stem = b_root

    h_leaf = FT(0.01) #10mm, guess
    h_stem = FT(18.5)# height of trunk, from Yujie's paper


    Kmax = 1.8e-10 #m^3/m^2/s/Pa, from Natan (10 mol/s/m^2/MPa) 


    K_max_stem = FT(Kmax)
    K_max_root = FT(Kmax)

    SAI = FT(0.00242) # Basal area per ground area
    LAI = FT(4.2) # from Yujie's paper
    f_root_to_shoot = FT(1.0 / 5.0) # guess
    RAI = SAI * f_root_to_shoot # following CLM
    # currently hardcoded to match the soil coordinates. this has to
    # be fixed eventually.
    z_root_depths = reverse(-Array(1:1:10.0) ./ 10.0 * 2.0 .+ 0.2 / 2.0)
    z_bottom_stem = FT(0.0)# this is OK
    z_leaf = h_stem

    roots_domain = RootDomain{FT}(z_root_depths, [z_bottom_stem, z_leaf])
    function root_distribution(z::T) where {T}
        return T(1.0 / 0.5) * exp(z / T(0.5))
    end


    roots_ps = Roots.RootsParameters{FT}(
        a_root,
        b_root,
        a_stem,
        b_stem,
        h_stem,
        h_leaf,
        K_max_root,
        K_max_stem,
        LAI,
        RAI,
        SAI,
        root_distribution, # exponential root distribution
    )

    zmin = FT(-2.0)
    zmax = FT(0.0)
    nelements = 20
    soil_domain = Column(FT, zlim = (zmin, zmax), nelements = nelements)
    ν = FT(0.495)
    Ksat = FT(4e-7) # matches Natan, m/s
    S_s = FT(1e-3) #inverse meters, guess
    vg_n = FT(1.5)
    vg_α = FT(0.10682) # inverse meters. From Natan (10.9/MPa)
    vg_m = FT(1) - FT(1) / vg_n
    θ_r = FT(0.0)
    soil_ps = Soil.RichardsParameters{FT}(ν, vg_α, vg_n, vg_m, Ksat, S_s, θ_r)

    soil_args = (
        domain = soil_domain,
        param_set = soil_ps,
        boundary_conditions = FluxBC(0.0, 0.0),
    )
    root_args = (domain = roots_domain, param_set = roots_ps)
    land_args = (precipitation = (t) -> 0.0, transpiration = (t) -> 0.0)

    land = RootSoilModel{FT}(;
        land_args = land_args,
        soil_model_type = Soil.RichardsModel{FT},
        soil_args = soil_args,
        vegetation_model_type = Roots.RootsModel{FT},
        vegetation_args = root_args,
    )
    Y, p, cds = initialize(land)
    ode! = make_ode_function(land)
    p_stem_ini = -0.5e6
    p_leaf_ini = -1e6
    θ_stem_0 = Roots.p_to_θ(p_stem_ini)
    θ_leaf_0 = Roots.p_to_θ(p_leaf_ini)
    Y.vegetation.θ .= FT.([θ_stem_0, θ_leaf_0])
    Y.soil.ϑ_l .= FT(0.4)
    update_aux! = make_update_aux(land)
    update_aux!(p, Y, 0.0)

    #sim
    t0 = FT(0)
    N_days = 1
    tf = FT(3600 * 24 * N_days)
    dt = FT(1)

    sv = SavedValues(FT, ClimaCore.Fields.FieldVector)
    daily = Array(2:(3600 * 24):(N_days * 3600 * 24))
    cb = SavingCallback(
        (u, t, integrator) -> copy(integrator.p),
        sv;
        saveat = daily,
    )
    prob = ODEProblem(ode!, Y, (t0, tf), p)
    sol = solve(prob, RK4(), dt = dt, callback = cb)
    #Currently just testing to make sure it runs, but need to have a better test suite.
end

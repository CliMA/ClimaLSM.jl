using Test
using UnPack
using DifferentialEquations
using OrdinaryDiffEq: ODEProblem, solve, Euler
using ClimaCore

if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
using ClimaLSM
using ClimaLSM.Domains: LSMMultiColumnDomain
using ClimaLSM.Soil
using ClimaLSM.Pond

FT = Float64
#@testset "Pond soil LSM integretation test" begin

    function precipitation(t::ft) where {ft}
        if t < ft(20)
            precip = -ft(1e-8)
        else
            precip = t < ft(100) ? -ft(5e-5) : ft(0.0)
        end
        return precip
    end
    ν = FT(0.495)
    Ksat = FT(0.0443 / 3600 / 100) # m/s
    S_s = FT(1e-3) #inverse meters
    vg_n = FT(2.0)
    vg_α = FT(2.6) # inverse meters
    vg_m = FT(1) - FT(1) / vg_n
    θ_r = FT(0)
    xmax = FT(10)
    xmin = FT(0)
    ymax = FT(10)
    ymin = FT(0)
    zmax = FT(0)
    zmin = FT(-1)
    nelems = (2,2,20)
    npoly = 1
    lsm_domain = LSMMultiColumnDomain(;xlim = (xmin, xmax), ylim = (ymin, ymax), zlim = (zmin, zmax), nelements = nelems, npolynomial = npoly)

    soil_domain = lsm_domain.subsurface
    soil_ps = Soil.RichardsParameters{FT}(ν, vg_α, vg_n, vg_m, Ksat, S_s, θ_r)
    soil_args = (domain = soil_domain, param_set = soil_ps)
    surface_water_args = (domain = lsm_domain.surface,)

    land_args = (precip = precipitation,)

    land = LandHydrology{FT}(;
        land_args = land_args,
        soil_model_type = Soil.RichardsModel{FT},
        soil_args = soil_args,
        surface_water_model_type = Pond.PondModel{FT},
        surface_water_args = surface_water_args,
    )
    Y, p, coords = initialize(land)
    function init_soil!(Ysoil, coords, params)
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
        Ysoil.soil.ϑ_l .= hydrostatic_profile.(coords.z, Ref(params))
    end
    init_soil!(Y, coords.soil, land.soil.param_set)
    # initialize the pond height to zero
    Y.surface_water.η .= 0.0

    saved_values = SavedValues(FT, ClimaCore.Fields.FieldVector)
    ode! = make_ode_function(land)
    t0 = FT(0)
    tf = FT(200)
    dt = FT(1)
    cb = SavingCallback((u, t, integrator) -> integrator.p, saved_values)
    prob = ODEProblem(ode!, Y, (t0, tf), p)
    sol = solve(prob, Euler(), dt = dt, callback = cb)
    # it running is the test

    # Infiltration at point test

    η = [0.0, 0.0, 0.0, 0.0, 1.0, 1.0]
    i_c = [2.0, 2.0, -2.0, -2.0, 2.0, 2.0]
    P = [3.0, -1.0, 3.0, -3.0, 1.0, 3.0]
    iap = [2.0, -1.0, -2.0, -3.0, 2.0, 2.0]
    @test infiltration_at_point.(η, i_c, P) ≈ iap
#end
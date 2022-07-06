using Test

using OrdinaryDiffEq: ODEProblem, solve, RK4, Euler
using DiffEqCallbacks
using DiffEqBase
using ClimaCore
import CLIMAParameters as CP

if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
using ClimaLSM.Bucket:
    BucketModel,
    BucketModelParameters,
    PrescribedAtmosphere,
    PrescribedRadiativeFluxes,
    BulkAlbedo
using ClimaLSM.Domains: coordinates, Plane
using ClimaLSM: initialize, make_update_aux, make_ode_function


FT = Float64

# Bucket model parameters
import ClimaLSM
include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))
earth_param_set = create_lsm_parameters(FT)

α_soil = (coordinate_point) -> FT(0.2) # soil albedo, spatially constant
α_snow = FT(0.8) # snow albedo
albedo = BulkAlbedo{FT}(α_snow, α_soil)
S_c = FT(0.2)
W_f = FT(0.15)
d_soil = FT(100.0) # soil depth
T0 = FT(280.0)
z_0m = FT(1e-2)
z_0b = FT(1e-3)
κ_soil = FT(1.5)
ρc_soil = FT(2e6)
bucket_parameters = BucketModelParameters(
    d_soil,
    T0,
    κ_soil,
    ρc_soil,
    albedo,
    S_c,
    W_f,
    z_0m,
    z_0b,
    earth_param_set,
)

# Model domain
bucket_domain = Plane(;
    xlim = (0.0, 1.0),
    ylim = (0.0, 1.0),
    nelements = (1, 1),
    periodic = (true, true),
    npolynomial = 1,
)

@testset "Zero flux RHS" begin
    "Radiation"
    SW_d = (t) -> eltype(t)(0.0)
    LW_d = (t) -> eltype(t)(5.67e-8 * 280.0^4.0)
    bucket_rad = PrescribedRadiativeFluxes(FT, SW_d, LW_d)
    "Atmos"
    precip = (t) -> eltype(t)(0) # no precipitation
    T_atmos = (t) -> eltype(t)(280.0)
    u_atmos = (t) -> eltype(t)(1.0)
    q_atmos = (t) -> eltype(t)(0.0) # no atmos water
    h_atmos = FT(30)
    ρ_atmos = (t) -> eltype(t)(1.13)
    ρ_sfc = FT(1.15)
    bucket_atmos = PrescribedAtmosphere(
        precip,
        T_atmos,
        u_atmos,
        q_atmos,
        ρ_atmos,
        h_atmos,
        ρ_sfc,
    )

    model = BucketModel(
        parameters = bucket_parameters,
        domain = bucket_domain,
        atmosphere = bucket_atmos,
        radiation = bucket_rad,
    )
    # Initial conditions with no moisture
    Y, p, coords = initialize(model)
    Y.bucket.T_sfc .= 280.0
    Y.bucket.W .= 0.0 # no moisture
    Y.bucket.Ws .= 0.0 # no runoff
    Y.bucket.S .= 0.0

    ode_function! = make_ode_function(model)
    dY = similar(Y)
    # init to nonzero numbers
    dY.bucket.T_sfc .= 1.0
    dY.bucket.W .= 1.0
    dY.bucket.Ws .= 1.0
    ode_function!(dY, Y, p, 0.0)
    @test sum(parent(dY.bucket.T_sfc)) < eps(FT)
    @test sum(parent(dY.bucket.W)) < eps(FT)
    @test sum(parent(dY.bucket.Ws)) < eps(FT)
end

@testset "Energy Only Equilibrium" begin
    # Energy: P = E = W = 0 (no moisture)
    # Set u_atmos = 0.0 so there are no turbulent fluxes
    # Set downward radiation = 300 W / m2
    # RHS of heat equation is -1/d(R↓-σT^4 +κ/d *(T-T_0))
    # can solve for fixed point = 255.69458073555182
    "Radiation"
    SW_d = (t) -> eltype(t)(290.0)
    LW_d = (t) -> eltype(t)(10.0)
    bucket_rad = PrescribedRadiativeFluxes(FT, SW_d, LW_d)
    # Atmosphere with no precipitation
    precip = (t) -> eltype(t)(0) # no precipitation
    T_atmos = (t) -> eltype(t)(298.0)
    u_atmos = (t) -> eltype(t)(0.0) # no turbulent fluxes
    q_atmos = (t) -> eltype(t)(0.0) # no atmos water
    h_atmos = FT(30)
    ρ_atmos = (t) -> eltype(t)(1.13)
    ρ_sfc = FT(1.15)
    bucket_atmos = PrescribedAtmosphere(
        precip,
        T_atmos,
        u_atmos,
        q_atmos,
        ρ_atmos,
        h_atmos,
        ρ_sfc,
    )

    model = BucketModel(
        parameters = bucket_parameters,
        domain = bucket_domain,
        atmosphere = bucket_atmos,
        radiation = bucket_rad,
    )

    # Initial conditions with no moisture
    Y, p, coords = initialize(model)
    Y.bucket.T_sfc .= 255.69458073555182
    Y.bucket.W .= 0.0 # no moisture
    Y.bucket.Ws .= 0.0 # no runoff
    Y.bucket.S .= 0.0
    ode_function! = make_ode_function(model)
    dY = similar(Y)
    ode_function!(dY, Y, p, 0.0)
    @test sum(parent(dY.bucket.T_sfc)) < eps(FT)
    @test sum(parent(dY.bucket.W)) < eps(FT)
    @test sum(parent(dY.bucket.Ws)) < eps(FT)
    @test sum(
        parent(
            p.bucket.R_n .-
            -(290.0 * 0.8 + 10.0 - 5.67e-8 * 255.69458073555182^4.0),
        ),
    ) < eps(FT)


end

@testset "Moisture Conservation" begin
    # Ensure no heat fluxes (no evaporation too, so no LHF)
    # and check that runoff+soil water balance precip.
    "Radiation"
    SW_d = (t) -> eltype(t)(0.0)
    LW_d = (t) -> eltype(t)(280.0^4.0 * 5.67e-8) # zeros out LW radiation
    bucket_rad = PrescribedRadiativeFluxes(FT, SW_d, LW_d)
    # Atmosphere with no precipitation
    precip = (t) -> eltype(t)(1e-2)
    T_atmos = (t) -> eltype(t)(298.0)
    u_atmos = (t) -> eltype(t)(0.0) # no turbulent fluxes
    q_atmos = (t) -> eltype(t)(0.0) # no atmos water
    h_atmos = FT(30)
    ρ_atmos = (t) -> eltype(t)(1.13)
    ρ_sfc = FT(1.15)
    bucket_atmos = PrescribedAtmosphere(
        precip,
        T_atmos,
        u_atmos,
        q_atmos,
        ρ_atmos,
        h_atmos,
        ρ_sfc,
    )

    model = BucketModel(
        parameters = bucket_parameters,
        domain = bucket_domain,
        atmosphere = bucket_atmos,
        radiation = bucket_rad,
    )

    # Initial conditions with no moisture
    Y, p, coords = initialize(model)
    Y.bucket.T_sfc .= 280.0
    Y.bucket.W .= 0.14 # no moisture
    Y.bucket.Ws .= 0.0 # no runoff
    Y.bucket.S .= 0.0
    ode_function! = make_ode_function(model)

    t0 = 0.0
    tf = 10.0
    prob = ODEProblem(ode_function!, Y, (t0, tf), p)
    Δt = 1.0
    integrator = init(prob, RK4(); dt = Δt, saveat = 0:Δt:tf)
    for step in 1:10
        step!(integrator, Δt, true)
    end

    # net soil water
    land_water = unique(
        parent(
            integrator.sol.u[end].bucket.W .+ integrator.sol.u[end].bucket.Ws,
        ),
    )[1]
    dland_water = land_water - 0.14 # 0.14 is initial water
    datmos_water = 1e-2 * integrator.sol.t[end]
    @test datmos_water ≈ dland_water
end

@testset "Energy + Moisture Conservation" begin
    "Radiation"
    SW_d = (t) -> eltype(t)(10.0)
    LW_d = (t) -> eltype(t)(300.0)
    bucket_rad = PrescribedRadiativeFluxes(FT, SW_d, LW_d)
    "Atmos"
    precip = (t) -> eltype(t)(1e-6)
    T_atmos = (t) -> eltype(t)(298.0)
    u_atmos = (t) -> eltype(t)(4.0)
    q_atmos = (t) -> eltype(t)(0.01)
    h_atmos = FT(30)
    ρ_atmos = (t) -> eltype(t)(1.13)
    ρ_sfc = FT(1.15)
    bucket_atmos = PrescribedAtmosphere(
        precip,
        T_atmos,
        u_atmos,
        q_atmos,
        ρ_atmos,
        h_atmos,
        ρ_sfc,
    )

    model = BucketModel(
        parameters = bucket_parameters,
        domain = bucket_domain,
        atmosphere = bucket_atmos,
        radiation = bucket_rad,
    )

    # Initial conditions with no moisture
    Y, p, coords = initialize(model)
    Y.bucket.T_sfc .= 280.0
    Y.bucket.W .= 0.149
    Y.bucket.Ws .= 0.0
    Y.bucket.S .= 0.0
    ode_function! = make_ode_function(model)

    t0 = 0.0
    tf = 10000.0
    prob = ODEProblem(ode_function!, Y, (t0, tf), p)
    Δt = 100.0
    # need a callback to get and store `p`
    saved_values = SavedValues(FT, ClimaCore.Fields.FieldVector)
    cb = SavingCallback(
        (u, t, integrator) -> copy(integrator.p),
        saved_values;
        saveat = 0:Δt:tf,
    )
    sol = solve(prob, Euler(); dt = Δt, saveat = 0:Δt:tf, callback = cb)
    # net soil water
    W = [unique(parent(sol.u[k].bucket.W))[1] for k in 1:length(sol.t)]
    Ws = [unique(parent(sol.u[k].bucket.Ws))[1] for k in 1:length(sol.t)]
    T_sfc = [unique(parent(sol.u[k].bucket.T_sfc))[1] for k in 1:length(sol.t)]
    R_n = [
        unique(parent(saved_values.saveval[k].bucket.R_n))[1] for
        k in 1:length(sol.t)
    ]
    SHF = [
        unique(parent(saved_values.saveval[k].bucket.SHF))[1] for
        k in 1:length(sol.t)
    ]
    LHF = [
        unique(parent(saved_values.saveval[k].bucket.LHF))[1] for
        k in 1:length(sol.t)
    ]
    E = [
        unique(parent(saved_values.saveval[k].bucket.E))[1] for
        k in 1:length(sol.t)
    ]

    F_g = -κ_soil .* (T_sfc .- 280.0) ./ d_soil
    F_sfc = LHF .+ SHF .+ R_n

    # dY(t+dt) = Y(t+dt) - W(t) = RHS|_t * dt for Euler stepping
    # First element of `p` not filled in, so start at 2.
    # If we call update aux before the simulation, it should be.


    # Look at fractional error
    dE_land = ρc_soil * (T_sfc[end] - T_sfc[2])
    dE_expected = -1 / d_soil * sum(F_sfc[2:(end - 1)] .- F_g[2:(end - 1)]) * Δt
    @test (dE_land - dE_expected) ./ (ρc_soil * sum(T_sfc) ./ length(sol.t)) <
          1e-14
    dW_land = W[end] + Ws[end] - W[2] - Ws[2]
    dW_expected = sum(precip.(sol.t[2:(end - 1)]) .- E[2:(end - 1)]) * Δt

    @test (dW_land - dW_expected) ./ 0.15 < 1e-14

end

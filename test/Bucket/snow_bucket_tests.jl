using Test

using OrdinaryDiffEq: ODEProblem, solve, Euler
using DiffEqCallbacks
using Statistics
using ClimaCore

if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
using ClimaLSM.Bucket:
    BucketModel,
    BucketModelParameters,
    PrescribedAtmosphere,
    PrescribedRadiativeFluxes,
    BulkAlbedo,
    partition_surface_fluxes
using ClimaLSM.Domains: coordinates, Plane
using ClimaLSM: initialize, make_update_aux, make_ode_function

FT = Float64

# Bucket model parameters
import ClimaLSM
import ClimaLSM.Parameters as LSMP
include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))
earth_param_set = create_lsm_parameters(FT)
α_soil = (coordinate_point) -> FT(0.2) # soil albedo, spatially constant
α_snow = FT(0.8) # snow albedo
albedo = BulkAlbedo{FT}(α_snow, α_soil)
σS_c = FT(0.2)
W_f = FT(0.15)
d_soil = FT(100.0) # soil depth
T0 = FT(280.0)
z_0m = FT(1e-2)
z_0b = FT(1e-3)
κ_soil = FT(1.5)
ρc_soil = FT(2e6)

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
    LW_d = (t) -> eltype(t)(5.67e-8 * 270.0^4.0)
    bucket_rad = PrescribedRadiativeFluxes(FT, SW_d, LW_d)
    "Atmos"
    precip = (t) -> eltype(t)(0) # no precipitation
    T_atmos = (t) -> eltype(t)(270.0)
    u_atmos = (t) -> eltype(t)(1.0)
    q_atmos = (t) -> eltype(t)(0.003280658023979249)
    h_atmos = FT(30)
    ρ_atmos = (t) -> eltype(t)(1.13)
    ρ_sfc = FT(1.15)
    bucket_atmos = PrescribedAtmosphere(
        precip,
        precip,
        T_atmos,
        u_atmos,
        q_atmos,
        ρ_atmos,
        h_atmos,
        ρ_sfc,
    )
    Δt = FT(1.0)
    bucket_parameters = BucketModelParameters(
        d_soil,
        T0,
        κ_soil,
        ρc_soil,
        albedo,
        σS_c,
        W_f,
        z_0m,
        z_0b,
        Δt,
        earth_param_set,
    )
    model = BucketModel(
        parameters = bucket_parameters,
        domain = bucket_domain,
        atmosphere = bucket_atmos,
        radiation = bucket_rad,
    )
    # Initial conditions with no moisture
    Y, p, coords = initialize(model)
    Y.bucket.T_sfc .= 270.0
    Y.bucket.W .= 0.0 # no moisture
    Y.bucket.Ws .= 0.0 # no runoff
    Y.bucket.σS .= 0.5

    ode_function! = make_ode_function(model)
    dY = similar(Y)
    # init to nonzero numbers
    dY.bucket.T_sfc .= 1.0
    dY.bucket.W .= 1.0
    dY.bucket.Ws .= 1.0
    dY.bucket.σS .= 1.0
    ode_function!(dY, Y, p, 0.0)
    @test abs(
        unique(parent(dY.bucket.T_sfc))[1] +
        (270.0 - 280.0) / 100.0^2.0 / ρc_soil * κ_soil,
    ) < eps(FT)
    @test sum(parent(dY.bucket.W)) < eps(FT)
    @test sum(parent(dY.bucket.Ws)) < eps(FT)
    @test sum(parent(dY.bucket.σS)) < eps(FT)
end

@testset "Melting only" begin
    "Radiation"
    SW_d = (t) -> eltype(t)(20.0)
    LW_d = (t) -> eltype(t)(5.67e-8 * 274.0^4.0)
    bucket_rad = PrescribedRadiativeFluxes(FT, SW_d, LW_d)
    "Atmos"
    precip = (t) -> eltype(t)(0) # no precipitation
    T_atmos = (t) -> eltype(t)(274.0)
    u_atmos = (t) -> eltype(t)(1.0)
    q_atmos = (t) -> eltype(t)(0.004506085305835068)
    h_atmos = FT(30)
    ρ_atmos = (t) -> eltype(t)(1.13)
    ρ_sfc = FT(1.15)
    bucket_atmos = PrescribedAtmosphere(
        precip,
        precip,
        T_atmos,
        u_atmos,
        q_atmos,
        ρ_atmos,
        h_atmos,
        ρ_sfc,
    )
    Δt = FT(1.0)
    bucket_parameters = BucketModelParameters(
        d_soil,
        T0,
        κ_soil,
        ρc_soil,
        albedo,
        σS_c,
        W_f,
        z_0m,
        z_0b,
        3.0 * Δt,
        earth_param_set,
    )
    model = BucketModel(
        parameters = bucket_parameters,
        domain = bucket_domain,
        atmosphere = bucket_atmos,
        radiation = bucket_rad,
    )
    # Initial conditions with no moisture
    Y, p, coords = initialize(model)
    Y.bucket.T_sfc .= 274.0
    Y.bucket.W .= 0.0 # no moisture
    Y.bucket.Ws .= 0.0 # no runoff
    Y.bucket.σS .= 0.5

    ode_function! = make_ode_function(model)
    dY = similar(Y)
    # init to nonzero numbers
    dY.bucket.T_sfc .= 1.0
    dY.bucket.W .= 1.0
    dY.bucket.Ws .= 1.0
    dY.bucket.σS .= 1.0
    ode_function!(dY, Y, p, 0.0)
    ΔW = unique(
        parent(dY.bucket.W .+ dY.bucket.σS .+ dY.bucket.Ws .+ p.bucket.E),
    )[1]
    _LH_f0 = LSMP.LH_f0(earth_param_set)
    # dI/dt = -F_snow
    ∑F_snow = unique(parent(1.0 .* dY.bucket.σS .* _LH_f0))[1]
    F_sfc = unique(parent(p.bucket.R_n .+ p.bucket.SHF .+ p.bucket.LHF))[1]
    G = F_sfc - ∑F_snow
    ΔE_soil = unique(parent(dY.bucket.T_sfc))[1] * ρc_soil
    F_bot_soil =
        unique(parent(@. -κ_soil / (d_soil) * (Y.bucket.T_sfc - T0)))[1]
    @test abs(ΔW) < eps(FT)
    @test abs(-(G .- F_bot_soil) / d_soil - ΔE_soil) < eps(FT)
    @test abs(F_sfc - ∑F_snow) < 10.0 * eps(FT)
end

@testset "Conservation of water and energy" begin
    "Radiation"
    SW_d = (t) -> eltype(t)(20.0)
    LW_d = (t) -> eltype(t)(20.0 + 5.67e-8 * 274.0^4.0)
    bucket_rad = PrescribedRadiativeFluxes(FT, SW_d, LW_d)
    "Atmos"
    precip = (t) -> eltype(t)(0) # no precipitation
    T_atmos = (t) -> eltype(t)(274.0)
    u_atmos = (t) -> eltype(t)(10.0)
    q_atmos = (t) -> eltype(t)(0.0)
    h_atmos = FT(30)
    ρ_atmos = (t) -> eltype(t)(1.13)
    ρ_sfc = FT(1.15)
    bucket_atmos = PrescribedAtmosphere(
        precip,
        precip,
        T_atmos,
        u_atmos,
        q_atmos,
        ρ_atmos,
        h_atmos,
        ρ_sfc,
    )
    Δt = FT(1.0)
    τc = FT(10.0)
    bucket_parameters = BucketModelParameters(
        d_soil,
        T0,
        κ_soil,
        ρc_soil,
        albedo,
        σS_c,
        W_f,
        z_0m,
        z_0b,
        τc,
        earth_param_set,
    )
    model = BucketModel(
        parameters = bucket_parameters,
        domain = bucket_domain,
        atmosphere = bucket_atmos,
        radiation = bucket_rad,
    )

    # run for some time
    Y, p, coords = initialize(model)
    Y.bucket.T_sfc .= 274.0
    Y.bucket.W .= 0.0 # no moisture
    Y.bucket.Ws .= 0.0 # no runoff
    Y.bucket.σS .= 0.5

    ode_function! = make_ode_function(model)
    update_aux! = make_update_aux(model)
    update_aux!(p, Y, 0.0)
    saved_values = SavedValues(FT, ClimaCore.Fields.FieldVector)
    cb = SavingCallback(
        (u, t, integrator) -> copy(integrator.p),
        saved_values;
        saveat = 0:Δt:100.0,
    )
    prob = ODEProblem(ode_function!, Y, (0.0, 100.0), p)
    sol = solve(
        prob,
        Euler();
        dt = Δt,
        reltol = 1e-6,
        abstol = 1e-6,
        callback = cb,
    )

    W = [unique(parent(sol.u[k].bucket.W))[1] for k in 1:length(sol.t)]
    Ws = [unique(parent(sol.u[k].bucket.Ws))[1] for k in 1:length(sol.t)]
    σS = [unique(parent(sol.u[k].bucket.σS))[1] for k in 1:length(sol.t)]
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
    _LH_f0 = LSMP.LH_f0(model.parameters.earth_param_set)
    _T_freeze = LSMP.T_freeze(model.parameters.earth_param_set)
    Isnow = -_LH_f0 * σS
    dIsnowdt = (Isnow[2:end] - Isnow[1:(end - 1)]) / Δt
    snow_cover_fraction(σS) = σS > eps(FT) ? FT(1.0) : FT(0.0)
    scf = snow_cover_fraction.(σS)
    partitioned_fluxes =
        partition_surface_fluxes.(
            σS,
            T_sfc,
            model.parameters.τc,
            scf,
            E,
            F_sfc,
            _LH_f0,
            _T_freeze,
        )
    F_melt =
        [partitioned_fluxes[k].F_melt for k in 1:length(partitioned_fluxes)]
    F_into_snow = [
        partitioned_fluxes[k].F_into_snow for k in 1:length(partitioned_fluxes)
    ]
    G = [partitioned_fluxes[k].G for k in 1:length(partitioned_fluxes)]

    # First point is off (why!), but the rest are ≈ 0

    # dIsnow/dt = -F_into_snow
    Δsnow = dIsnowdt .+ F_into_snow[1:(end - 1)]
    @test maximum(abs.(Δsnow[2:end])) ./ 0.5 / _LH_f0 * Δt < eps(FT)
    # dT/dt = -1/c*1/d*(G-F_g)
    dTdt = (T_sfc[2:end] - T_sfc[1:(end - 1)]) / Δt
    ΔT = dTdt .* ρc_soil * d_soil .+ (G[1:(end - 1)] .- F_g[1:(end - 1)])
    @test maximum(abs.(ΔT[2:end])) ./ (274.0 * ρc_soil * d_soil) * Δt < eps(FT)
    # dW/dt = -E -> dWdt + E = 0
    WL = W .+ σS .+ Ws
    dWdt = (WL[2:end] - WL[1:(end - 1)]) / Δt

    @test maximum(abs.(dWdt[2:end] .+ E[2:(end - 1)])) < eps(FT)

    IL = -_LH_f0 * σS .+ d_soil * ρc_soil .* T_sfc
    dILdt = (IL[2:end] - IL[1:(end - 1)]) / Δt
    # dIL/dt = - (F_sfc - F_g)
    Δ = dILdt .+ (F_sfc[1:(end - 1)] .- F_g[1:(end - 1)])
    @test maximum(abs.(Δ[2:end])) ./ (274.0 * ρc_soil * d_soil) * Δt <
          2.0 * eps(FT)

end

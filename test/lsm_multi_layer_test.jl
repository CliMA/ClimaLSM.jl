using Test
using UnPack
using NLsolve
using OrdinaryDiffEq: ODEProblem, solve, SSPRK33
using DifferentialEquations
using ClimaCore
if !("." in LOAD_PATH) # for ease of include
    push!(LOAD_PATH, ".")
end

using ClimaLSM
using ClimaLSM.Domains: Column, RootDomain
using ClimaLSM.Soil
using ClimaLSM.Roots



FT = Float64
saved_values = SavedValues(FT, ClimaCore.Fields.FieldVector)

const a_root = FT(13192)
const a_stem = FT(515.5605)
const b_root = FT(2.1079/1e6) # Inverse Pa
const b_stem = FT(0.9631/1e6) # Inverse Pa
const h_leaf = FT(0.01) #10mm, guess
const h_stem = FT(0.5)
const K_max_stem = FT(3.75e-9)
const K_max_root = FT(1.48e-8)
const SAI = FT(0.1)
const RAI = FT(0.1)
const LAI = FT(0.3)
# currently hardcoded to match the soil coordinates. this has to
# be fixed eventually.
const z_root_depths = -Array(1:1:20.0) ./ 20.0 * 3.0 .+ 0.15 / 2.0
const z_bottom_stem = FT(0.0)
const z_leaf = FT(0.5) # height of leaf

roots_domain = RootDomain{FT}(z_root_depths, [z_bottom_stem, z_leaf])

function root_distribution(z::T) where {T}
    return  T(1.0/0.25)*exp(z/T(0.25))
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
    root_distribution # exponential root distribution
)

zmin = FT(-3.0)
zmax = FT(0.0)
nelements = 20
soil_domain = Column(FT, zlim = (zmin, zmax), nelements = nelements);
const ν = FT(0.495);
const Ksat = FT(0.0443 / 3600 / 100); # m/s
const S_s = FT(1e-3); #inverse meters
const vg_n = FT(2.0);
const vg_α = FT(2.6); # inverse meters
const vg_m = FT(1) - FT(1) / vg_n;
const θ_r = FT(0);
soil_ps = Soil.RichardsParameters{FT}(ν, vg_α, vg_n, vg_m, Ksat, S_s, θ_r);

soil_args = (domain = soil_domain, param_set = soil_ps)
root_args = (domain = roots_domain, param_set = roots_ps)
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
init_soil!(Y, coords.soil, land.soil.param_set)

## soil is at total ψ+z = -3.0 #m
## Want ρgΨ_plant = ρg(-3) - ρg z_plant
# we should standardize the units! and not ahve to convert every time.
# convert parameters once up front and then not each RHS
p_stem_ini = -28665.0
p_leaf_ini = 0.0

theta_stem_0 = p_to_theta(p_stem_ini)
theta_leaf_0 = p_to_theta(p_leaf_ini)
y0 = FT.([theta_stem_0, theta_leaf_0])
Y.vegetation.theta .= y0

ode! = make_ode_function(land)
t0 = FT(0);
tf = FT(2)
dt = FT(1);
cb = SavingCallback((u, t, integrator) -> integrator.p, saved_values)
prob = ODEProblem(ode!, Y, (t0, tf), p);
sol = solve(prob, Euler(), dt = dt, callback = cb);
#Currently just testing to make sure it runs, but need to have a better test suite.
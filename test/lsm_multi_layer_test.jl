using Test
using UnPack
using NLsolve
using OrdinaryDiffEq: ODEProblem, solve, Euler
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
const b_root = FT(2.1079)
const b_stem = FT(0.9631)
const size_reservoir_leaf_moles = FT(16766.2790)
const size_reservoir_stem_moles = FT(11000.8837)
const K_max_root_moles = FT(12.9216)
const K_max_stem_moles = FT(3.4415)
const z_leaf = FT(12) # height of leaf
# currently hardcoded to match the soil coordinates. this has to
# be fixed eventually.
const z_root_depths = -Array(1:1:20.0) ./ 20.0 * 3.0 .+ 0.15 / 2.0
const z_bottom_stem = FT(0.0)
roots_domain = RootDomain{FT}(z_root_depths, [z_bottom_stem, z_leaf])
roots_ps = Roots.RootsParameters{FT}(
    a_root,
    b_root,
    a_stem,
    b_stem,
    size_reservoir_stem_moles,
    size_reservoir_leaf_moles,
    K_max_root_moles,
    K_max_stem_moles,
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
t0 = FT(0);
tf = FT(2)
dt = FT(1);
cb = SavingCallback((u, t, integrator) -> integrator.p, saved_values)
prob = ODEProblem(ode!, Y, (t0, tf), p);
sol = solve(prob, Euler(), dt = dt, callback = cb);
#Currently just testing to make sure it runs, but need to have a better test suite.

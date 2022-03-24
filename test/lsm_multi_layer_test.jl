using Test
using UnPack
using NLsolve
using OrdinaryDiffEq: ODEProblem, solve, Euler, Midpoint
using DifferentialEquations
using ClimaCore
if !("." in LOAD_PATH) # for ease of include
    push!(LOAD_PATH, ".")
end
using DelimitedFiles
using Dierckx
using Plots

using ClimaLSM
using ClimaLSM.Domains: Column, RootDomain
using ClimaLSM.Soil
using ClimaLSM.Roots

precip_θ_T = readdlm("/Users/katherinedeck/Desktop/ozark_site/p_data_2017.csv",',')
et = readdlm("/Users/katherinedeck/Desktop/ozark_site/et_data_2017.csv",',')
et_spline = Spline1D(et[:,1], et[:,2])
p_spline = Spline1D(precip_θ_T[:,1], precip_θ_T[:,2])
t = precip_θ_T[:,1]

precip_function(t::ft) where {ft} = p_spline(t)
transpiration_function(t::ft) where {ft} = et_spline(t)


FT = Float64

# Somewhat close to WEibull with C = 0.953 and B = 5.703 MPa
const a_root = FT(0.1)
const a_stem = a_root
const b_root = FT(0.17/1e6) # Inverse Pa
const b_stem = b_root

const h_leaf = FT(0.01) #10mm, guess
const h_stem = FT(18.5)# height of trunk, from Yujie's paper


Kmax =  7e-11*20.0 #m^3*m/m^2/s/Pa relative to BASAL area # conductance
# multiply by height to get conductivity

const K_max_stem = FT(Kmax)# ratio of 2:1:1 for roots:stem:leaves
const K_max_root = FT(Kmax)

const SAI = FT(0.00242) # Basal area per ground area
const LAI = FT(4.2) # from Yujie's paper
const f_root_to_shoot = FT(1.0/5.0) # guess
const RAI = SAI*f_root_to_shoot # following CLM
# currently hardcoded to match the soil coordinates. this has to
# be fixed eventually.
const z_root_depths = -Array(1:1:10.0) ./ 10.0 * 2.0 .+ 0.2 / 2.0# OK
const z_bottom_stem = FT(0.0)# this is OK
const z_leaf = h_stem 

roots_domain = RootDomain{FT}(z_root_depths, [z_bottom_stem, z_leaf])
# 0.95 is from CLM
function root_distribution(z::T) where {T}
    if z > -0.1
        return 1.0
    else
        return 0.0
    end
    
#    return  T(1.0/0.95)*exp(z/T(0.95))
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

zmin = FT(-2.0)
zmax = FT(0.0)
nelements = 10
soil_domain = Column(FT, zlim = (zmin, zmax), nelements = nelements);
const ν = FT(0.45);
const Ksat = FT(0.45 / 3600 / 100); # m/s,  0.45 cm/h
const S_s = FT(1e-3); #inverse meters, guess
const vg_n = FT(1.41);
const vg_α = FT(2.0); # inverse meters
const vg_m = FT(1) - FT(1) / vg_n;
const θ_r = FT(0.067);
soil_ps = Soil.RichardsParameters{FT}(ν, vg_α, vg_n, vg_m, Ksat, S_s, θ_r);

soil_args = (domain = soil_domain, param_set = soil_ps)
root_args = (domain = roots_domain, param_set = roots_ps)
land_args = (precipitation = 0.0, transpiration = (t) -> 0.0)

land = RootSoilModel{FT}(;
                         land_args = land_args,
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
        z_∇ = FT(-2)# matches zmin
        S = FT((FT(1) + (vg_α * (z - z_∇))^vg_n)^(-vg_m))
        ϑ_l = S * (ν - θ_r) + θ_r
        return FT(ϑ_l)
    end
    Ysoil.soil.ϑ_l .= hydrostatic_profile.(z, Ref(params))
end
init_soil!(Y, coords.soil, land.soil.param_set)
## soil is at total ψ+z = -2.0 #m
## Want ρgΨ_plant = ρg(-2) - ρg z_plant
# we should standardize the units! and not ahve to convert every time.
# convert parameters once up front and then not each RHS
p_stem_ini = -576715.0136522778
p_leaf_ini = -1e6

θ_stem_0 = Roots.p_to_θ(p_stem_ini)
θ_leaf_0 = Roots.p_to_θ(p_leaf_ini)
y0 = FT.([θ_stem_0, θ_leaf_0])
Y.vegetation.θ .= y0

ode! = make_ode_function(land)
t0 = FT(0);
tf = FT(3600*24*3.6)
dt = FT(100);
#update_aux! = make_update_aux(land)
#update_aux!(p,Y,t0)
sv = SavedValues(FT, ClimaCore.Fields.FieldVector)
cb = SavingCallback((u, t, integrator) -> copy(integrator.p), sv)
prob = ODEProblem(ode!, Y, (t0, tf), p);
#integrator = init(prob, Euler(); dt = dt)
sol = solve(prob, Euler(), dt = dt, callback = cb);
#Currently just testing to make sure it runs, but need to have a better test suite.

p_soil = [parent(sv.saveval[k].soil.ψ .+ coords.soil) .* 9800 for k in 1:1:length(sol.t)]

p_stem = [(Roots.θ_to_p.(sol.u[k].vegetation.θ[1]) .+ 0* 9800) for k in 1:1:length(sol.t)]
p_leaf = [(Roots.θ_to_p.(sol.u[k].vegetation.θ[2]) .+ 18.5.* 9800) for k in 1:1:length(sol.t)]
gaf = [sum(sv.saveval[k].root_extraction_source) for k in 1:1:length(sol.t)]

using Test
using Statistics
using DifferentialEquations
using UnPack
using NLsolve
using OrdinaryDiffEq: ODEProblem, solve, Euler
using ClimaCore
using DelimitedFiles
using Dierckx

if !("." in LOAD_PATH) # for ease of include
    push!(LOAD_PATH, ".")
end
using ClimaLSM
using ClimaLSM.Domains: Column, RootDomain
using ClimaLSM.Soil
using ClimaLSM.Roots

FT = Float64

const mass_mole_water = FT(0.0180) 
const a_root = FT(13192)
const a_stem = FT(515.5605)
const b_root = FT(2.1079) # Inverse MPa
const b_stem = FT(0.9631) # Inverse MPa
const h_stem = FT(13) # height of trunk
const K_max_stem = FT(3.4415*h_stem) 
const K_max_root = FT(601.1975)
const z_root_depths = Array(0:1:2) 
const z_bottom_stem = FT(5)
const z_leaf = h_stem 

root_domain = RootDomain{FT}(z_root_depths, [z_bottom_stem, z_leaf])

param_set = Roots.RootsParameters{FT}(
    a_root,
    b_root,
    a_stem,
    b_stem,
    h_stem,
    K_max_root,
    K_max_stem,
    )

function leaf_transpiration(t::ft) where {ft}
    T = ft(0.01/mass_mole_water) 
    return T
end

const p_soil0 = [-0.1; -0.2; -0.3] 
transpiration = PrescribedTranspiration{FT}((t::FT) -> leaf_transpiration(t))
root_extraction = PrescribedSoilPressure{FT}((t::FT) -> p_soil0)
roots = Roots.RootsModel{FT}(;
    domain = root_domain,
    param_set = param_set,
    root_extraction = root_extraction,
    transpiration = transpiration,
)
Y,p,_ = initialize(roots)
p_stem_ini = FT(-0.5)
p_leaf_ini = FT(-0.4) 

θ_stem_0 = Roots.p_to_θ(p_stem_ini)
θ_leaf_0 = Roots.p_to_θ(p_leaf_ini)
y0 = FT.([θ_stem_0, θ_leaf_0])
Y.vegetation.θ .= y0

ode! = make_ode_function(roots)
t0 = FT(0);
tf = FT(60);
dt = FT(1);
#update_aux! = make_update_aux(land)
#update_aux!(p,Y,t0)
sv = SavedValues(FT, ClimaCore.Fields.FieldVector)
cb = SavingCallback((u, t, integrator) -> copy(integrator.p), sv)
prob = ODEProblem(ode!, Y, (t0, tf), p);
#integrator = init(prob, Euler(); dt = dt)
sol = solve(prob, Euler(), dt = dt, callback = cb);
p_stem = [(Roots.θ_to_p.(sol.u[k].vegetation.θ[1])) for k in 1:1:length(sol.t)]
p_leaf = [(Roots.θ_to_p.(sol.u[k].vegetation.θ[2])) for k in 1:1:length(sol.t)]

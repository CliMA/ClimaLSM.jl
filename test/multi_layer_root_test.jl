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
const b_root = FT(2.1079) # Inverse Pa
const b_stem = FT(0.9631)

#const h_leaf = FT(0.01) #10mm, guess
const h_stem = FT(13)# height of trunk, from Yujie's paper

# Kmax =  7e-11*20.0 #m^3*m/m^2/s/Pa relative to BASAL area # conductance
# multiply by height to get conductivity

const K_max_stem = FT(3.4415*h_stem) # ratio of 2:1:1 for roots:stem:leaves
const K_max_root = FT(601.1975)

#const SAI = FT(0.00242) # Basal area per ground area
#const LAI = FT(4.2) # from Yujie's paper
#const f_root_to_shoot = FT(1.0/5.0) # guess
#const RAI = SAI*f_root_to_shoot # following CLM
# currently hardcoded to match the soil coordinates. this has to
# be fixed eventually.
const z_root_depths = Array(0:1:2) #reverse(-Array(1:1:10.0) ./ 10.0 * 2.0 .+ 0.2 / 2.0)
const z_bottom_stem = FT(5)# this is OK
const z_leaf = h_stem 

root_domain = RootDomain{FT}(z_root_depths, [z_bottom_stem, z_leaf])
#function root_distribution(z::T) where {T}
#    if z > -0.1
#        return 1.0
#    else
#        return 0.0
#    end
#    
#    return  T(1.0/0.95)*exp(z/T(0.95))
#end
param_set = Roots.RootsParameters{FT}(
    a_root,
    b_root,
    a_stem,
    b_stem,
    h_stem,
    K_max_root,
    K_max_stem,
    #LAI,
    #RAI,
    #SAI,
    #root_distribution, # exponential root distribution
)

function leaf_transpiration(t::ft) where {ft}
    T = ft(0.01/mass_mole_water) 
    return T
end

#const p_soil0 = (-2.0 .- z_root_depths) .* 9800
const p_soil0 = [-0.1; -0.2; -0.3] #zeros(10).-9800.0
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
tf = FT(60*60.0*2)
dt = FT(1);
#update_aux! = make_update_aux(land)
#update_aux!(p,Y,t0)
sv = SavedValues(FT, ClimaCore.Fields.FieldVector)
cb = SavingCallback((u, t, integrator) -> copy(integrator.p), sv)
prob = ODEProblem(ode!, Y, (t0, tf), p);
#integrator = init(prob, Euler(); dt = dt)
sol = solve(prob, Euler(), dt = dt, callback = cb);
p_stem = [(Roots.θ_to_p.(sol.u[k].vegetation.θ[1]) .+ 0* 9800) for k in 1:1:length(sol.t)]
p_leaf = [(Roots.θ_to_p.(sol.u[k].vegetation.θ[2]) .+ 18.5.* 9800) for k in 1:1:length(sol.t)]

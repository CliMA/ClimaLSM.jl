using Test
using Statistics
using DifferentialEquations
using UnPack
using NLsolve
using OrdinaryDiffEq: ODEProblem, solve, RK4
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
const a_root = FT(0.1)
const a_stem = a_root
const b_root = FT(0.17/1e6) # Inverse Pa
const b_stem = b_root

const h_leaf = FT(0.01) #10mm, guess
const h_stem = FT(18.5)# height of trunk, from Yujie's paper


Kmax =  1.8e-10 #m^3/m^2/s/Pa, from Natan (10 mol/s/m^2/MPa) 

const K_max_stem = FT(Kmax)
const K_max_root = FT(Kmax)

const SAI = FT(0.00242) # Basal area per ground area
const LAI = FT(4.2) # from Yujie's paper
const f_root_to_shoot = FT(1.0/5.0) # guess
const RAI = SAI*f_root_to_shoot # following CLM
# currently hardcoded to match the soil coordinates. this has to
# be fixed eventually.
const z_root_depths = reverse(-Array(1:1:10.0) ./ 10.0 * 2.0 .+ 0.2 / 2.0)
const z_bottom_stem = FT(0.0)# this is OK
const z_leaf = h_stem 

root_domain = RootDomain{FT}(z_root_depths, [z_bottom_stem, z_leaf])
function root_distribution(z::T) where {T}
    if z > -2.0
        return 1.0/2.0
    else
        return 0.0
    end
    
#    return  T(1.0/0.95)*exp(z/T(0.95))
end
param_set = Roots.RootsParameters{FT}(
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

function leaf_transpiration(t::ft) where {ft}
    T = ft(0.0)
 
    return T
end

#const p_soil0 = (-2.0 .- z_root_depths) .* 9800
const p_soil0 = zeros(10).-1e6
transpiration = PrescribedTranspiration{FT}((t::FT) -> leaf_transpiration(t))
root_extraction = PrescribedSoilPressure{FT}((t::FT) -> p_soil0)
roots = Roots.RootsModel{FT}(;
    domain = root_domain,
    param_set = param_set,
    root_extraction = root_extraction,
    transpiration = transpiration,
)
Y,p,_ = initialize(roots)
p_stem_ini = -1e6
p_leaf_ini = -3e6

θ_stem_0 = Roots.p_to_θ(p_stem_ini)
θ_leaf_0 = Roots.p_to_θ(p_leaf_ini)
y0 = FT.([θ_stem_0, θ_leaf_0])
Y.vegetation.θ .= y0

ode! = make_ode_function(roots)
t0 = FT(0);
tf = FT(3600*24*30)
dt = FT(1);
#update_aux! = make_update_aux(land)
#update_aux!(p,Y,t0)
sv = SavedValues(FT, ClimaCore.Fields.FieldVector)
cb = SavingCallback((u, t, integrator) -> copy(integrator.p), sv)
prob = ODEProblem(ode!, Y, (t0, tf), p);
#integrator = init(prob, Euler(); dt = dt)
sol = solve(prob, RK4(), dt = dt, callback = cb);
ϕ_stem = [(Roots.θ_to_p.(sol.u[k].vegetation.θ[1]) .+ 0* 9800) for k in 1:1:length(sol.t)];
ϕ_leaf = [(Roots.θ_to_p.(sol.u[k].vegetation.θ[2]) .+ 18.5.* 9800) for k in 1:1:length(sol.t)];

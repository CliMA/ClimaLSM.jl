using Test
using Statistics
using DifferentialEquations
using UnPack
using NLsolve
using OrdinaryDiffEq: ODEProblem, solve, Euler
using ClimaCore

if !("." in LOAD_PATH) # for ease of include
    push!(LOAD_PATH, ".")
end
using ClimaLSM
using ClimaLSM.Domains: Column, RootDomain
using ClimaLSM.Soil
using ClimaLSM.Roots


FT = Float64
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
const z_leaf = FT(0.5) # height of leaf
const z_root_depths = FT.([-0.4,-0.3,-0.2,-0.1]) # m, rooting depth
const z_bottom_stem = FT(0.0)

root_domain = RootDomain{FT}(z_root_depths, [z_bottom_stem, z_leaf])
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
    (z) -> FT(1.0/0.25)*exp(z/FT(0.25)) # exponential root distribution
)

function leaf_transpiration(t::ft) where {ft}
    T = ft(0.0)
    T_0 = FT(1e-5)
    if t < ft(500)
        T = T_0
    elseif t < ft(1000)
        T = ft(10 * (T_0 / 5) * (t - 500) / 500 + T_0)
    else
        T = ft(10 * (T_0 / 5) * 500 / 500 + T_0)
    end
    return T
end

const p_soil0 = FT.([-2e6, -2e6, -2e6, -2e6])
transpiration = PrescribedTranspiration{FT}((t::FT) -> leaf_transpiration(t))
root_extraction = PrescribedSoilPressure{FT}((t::FT) -> p_soil0)
roots = Roots.RootsModel{FT}(;
    domain = root_domain,
    param_set = param_set,
    root_extraction = root_extraction,
    transpiration = transpiration,
)

# Set system to equilibrium state by setting LHS of both ODEs to 0


function f!(F, P)
    T0 = 1e-5
    flow_in_stem =
        sum(
            ground_area_flux.(
                roots.domain.root_depths,
                roots.domain.compartment_heights[1],
                p_soil0,
                P[1],
                a_root,
                b_root,
                K_max_root,
                RAI,
            ) .* roots.param_set.root_distribution_function.(roots.domain.root_depths)
            .* (vcat(roots.domain.root_depths,[0.0])[2:end] - vcat(roots.domain.root_depths,[0.0])[1:end-1]))
    
    flow_out_stem = ground_area_flux(
        z_bottom_stem,
        z_leaf,
        P[1],
        P[2],
        a_stem,
        b_stem,
        K_max_stem,
        SAI,
    )
    F[1] = (flow_in_stem - flow_out_stem)
    F[2] = (flow_out_stem - T0*LAI)
end

soln = nlsolve(f!, [-1e6, -0.9e6], ftol = 1e-10)
p_stem_ini = soln.zero[1]
p_leaf_ini = soln.zero[2]

θ_stem_ini = p_to_θ(p_stem_ini)
θ_leaf_ini = p_to_θ(p_leaf_ini)
y0 = FT.([θ_stem_ini, θ_leaf_ini])

Y, p, coords = initialize(roots)
Y.vegetation.θ .= y0

root_ode! = make_ode_function(roots)

t0 = FT(0);
tf = FT(1200)
dt = FT(1);

prob = ODEProblem(root_ode!, Y, (t0, tf), p);
#integrator = init(prob, Euler(); dt = dt)
sol = solve(prob, Euler(), dt = dt);

dY = similar(Y)
root_ode!(dY, Y, p, 0.0)
@test sqrt(mean(dY.vegetation.θ .^ 2.0)) < 1e-8 # starts in equilibrium


y_1 = reduce(hcat, sol.u)[1, :]
y_2 = reduce(hcat, sol.u)[2, :]
p_stem = θ_to_p.(y_1)
p_leaf = θ_to_p.(y_2)

function f2!(F, P)
    p_soilf = p_soil0
    Tf = 1e-5* 3.0
    flow_in_stem =
        sum(
            ground_area_flux.(
                roots.domain.root_depths,
                roots.domain.compartment_heights[1],
                p_soilf,
                P[1],
                a_root,
                b_root,
                K_max_root,
                RAI,
            ) .* roots.param_set.root_distribution_function.(roots.domain.root_depths)
            .* (vcat(roots.domain.root_depths,[0.0])[2:end] - vcat(roots.domain.root_depths,[0.0])[1:end-1]))
    
    flow_out_stem = ground_area_flux(
        z_bottom_stem,
        z_leaf,
        P[1],
        P[2],
        a_stem,
        b_stem,
        K_max_stem,
        SAI,
    )
    F[1] = (flow_in_stem - flow_out_stem)
    F[2] = (flow_out_stem - Tf*LAI)

end


# Check that the final state is in the new equilibrium
soln = nlsolve(f2!, [p_leaf[end], p_stem[end]]; ftol = 1e-12)
p_stem_f = soln.zero[1]
p_leaf_f = soln.zero[2]
@test abs(p_stem_f - p_stem[end]) < 1e-7
@test abs(p_leaf_f - p_leaf[end]) < 1e-7

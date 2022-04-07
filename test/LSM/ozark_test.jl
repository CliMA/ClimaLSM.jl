using Test
using UnPack
using NLsolve
using OrdinaryDiffEq: ODEProblem, solve, Euler, RK4
using DifferentialEquations
using ClimaCore
if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
using DelimitedFiles
using Dierckx
using Plots
using Statistics
using Dates

using ClimaLSM
using ClimaLSM.Domains: Column, RootDomain
using ClimaLSM.Soil
using ClimaLSM.Roots

precip_θ_T =
    readdlm("/Users/katherinedeck/Desktop/ozark_site/p_data_2005.csv", ',')
et = readdlm("/Users/katherinedeck/Desktop/ozark_site/et_data_2005.csv", ',')
#et_spline = Spline1D(et[:, 1], max.(0.0, et[:, 2]))
p_spline = Spline1D(precip_θ_T[:, 1], -precip_θ_T[:, 2])
t = precip_θ_T[:, 1]



## Natan's data
data = readdlm(
    "/Users/katherinedeck/Downloads/holtzman_clima_output_april1.csv",
    ',',
)
natan_et = data[2:end, 15] * 18 / 1e3 / 1e3 # convert to m^3/m^2/s
swc_column = data[2:end, 20]
swc_surface = data[2:end, 19]
swc_obs_surface = data[2:end, 12]
lwp = data[2:end,18]
dates = data[2:end, 1]
data = nothing #frees memory
dates_julia = tryparse.(DateTime, dates)
our_year = dates_julia[Dates.year.(dates_julia) .== 2005]
seconds = Dates.value.(our_year .- our_year[1]) ./ 1000
our_year_swc_column = FT.(swc_column[Dates.year.(dates_julia) .== 2005])
our_year_swc_surface = FT.(swc_surface[Dates.year.(dates_julia) .== 2005])
our_year_swc_obs_surface = FT.(swc_obs_surface[Dates.year.(dates_julia) .== 2005])
our_year_lwp = FT.(lwp[Dates.year.(dates_julia) .== 2005])
our_year_et = FT.(natan_et[Dates.year.(dates_julia) .== 2005])
et_spline = Spline1D(seconds, our_year_et)



# Capping T >0 and P <0. Not sure. P>0 seemed to be fine, but ET <0 led to some issues with
# the integration (numerical instability)
precip_function(t::ft) where {ft} = p_spline(t) < 0.0 ? p_spline(t) : 0.0
transpiration_function(t::ft) where {ft} =
    et_spline(t) > 0.0 ? et_spline(t) : 0.0


FT = Float64

# Somewhat close to Weibull with C = 0.953 and B = 5.703 MPa
const a_root = FT(0.1)
const a_stem = a_root
const b_root = FT(0.17 / 1e6) # Inverse Pa
const b_stem = b_root
const h_leaf = FT(0.005) #5mm, guess
const h_stem = FT(18.5)# height of trunk, from Yujie's paper
Kmax = 1.8e-10 #m^3/m^2/s/Pa, from Natan (10 mol/s/m^2/MPa) 
const K_max_stem = FT(Kmax)
const K_max_root = FT(Kmax)
const SAI = FT(0.00242) # Basal area per ground area
const LAI = FT(4.2) # from Yujie's paper
const f_root_to_shoot = FT(1.0 / 5.0) # guess
const RAI = SAI * f_root_to_shoot # following CLM
# currently hardcoded to match the soil coordinates. this has to
# be fixed eventually.
const z_root_depths = reverse(-Array(1:1:10.0) ./ 10.0 * 2.0 .+ 0.2 / 2.0)
const z_bottom_stem = FT(0.0)# this is OK
const z_leaf = h_stem

roots_domain = RootDomain{FT}(z_root_depths, [z_bottom_stem, z_leaf])
# 0.5 is from Natan
function root_distribution(z::T) where {T}
    return T(1.0 / 0.5) * exp(z / T(0.5))
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
    root_distribution, # exponential root distribution
)

zmin = FT(-2.0)
zmax = FT(0.0)
nelements = 10
soil_domain = Column(FT, zlim = (zmin, zmax), nelements = nelements);
const ν = FT(0.55);
const Ksat = FT(4e-7) # matches Natan, m/s
const S_s = FT(1e-3); #inverse meters, guess
const vg_n = FT(1.5);
const vg_α = FT(0.10682); # inverse meters. From Natan (10.9/MPa)
const vg_m = FT(1) - FT(1) / vg_n;
const θ_r = FT(0.0);
soil_ps = Soil.RichardsParameters{FT}(ν, vg_α, vg_n, vg_m, Ksat, S_s, θ_r);

# To get the equilibrium IC
#=
soil_args = (domain = soil_domain, param_set = soil_ps, boundary_conditions = FluxBC(0.0,0.0))
root_args = (domain = roots_domain, param_set = roots_ps)

land_args = (precipitation = (t) -> 0.0, transpiration = (t) -> 0.0)

land = RootSoilModel{FT}(;
                         land_args = land_args,
                         soil_model_type = Soil.RichardsModel{FT},
                         soil_args = soil_args,
                         vegetation_model_type = Roots.RootsModel{FT},
                         vegetation_args = root_args,
                         )
Y, p, cds = initialize(land)
ode! = make_ode_function(land)
p_stem_ini = -0.5e5
p_leaf_ini = -1e5
θ_stem_0 = Roots.p_to_θ(p_stem_ini)
θ_leaf_0 = Roots.p_to_θ(p_leaf_ini)
Y.vegetation.θ .= FT.([θ_stem_0, θ_leaf_0])
Y.soil.ϑ_l .= FT(0.35)
update_aux! = make_update_aux(land)
update_aux!(p,Y,0.0)

#sim
t0 = FT(0);
N_days = 300
tf = FT(3600*24*N_days)
dt = FT(1);

sv = SavedValues(FT, ClimaCore.Fields.FieldVector)
daily = Array(2:3600*24:N_days*3600*24)
cb = SavingCallback((u, t, integrator) -> copy(integrator.p), sv; saveat = daily)
prob = ODEProblem(ode!, Y, (t0, tf), p);
sol = solve(prob, RK4(), dt = dt, callback = cb);

ϕ_soil =
    [parent(sv.saveval[k].soil.ψ .+ cds.soil.z) .* 9800 for k in 1:1:N_days]
pl = scatter()
function plot_each_day(pl, indices)
    for i in indices
        plot!(pl, ϕ_soil[i][:] ./ 1e6, parent(cds.soil.z), label = "")
    end
    plot!(
        pl,
        xlabel = "ϕ_soil (MPa)",
        ylabel = "z(m)",
        legend = :bottomright,
    )
end
plot_each_day(pl, 1:1:N_days)
=#
###

soil_args = (
    domain = soil_domain,
    param_set = soil_ps,
    boundary_conditions = PrecipFreeDrainage{FT}(precip_function),
)
root_args = (domain = roots_domain, param_set = roots_ps)
land_args =
    (precipitation = precip_function, transpiration = transpiration_function)

land = RootSoilModel{FT}(;
    land_args = land_args,
    soil_model_type = Soil.RichardsModel{FT},
    soil_args = soil_args,
    vegetation_model_type = Roots.RootsModel{FT},
    vegetation_args = root_args,
)
Y, p, cds = initialize(land)
ode! = make_ode_function(land)

# specify ICs from equilibrium run.
Y.vegetation.θ .=  [0.999596628794538, 0.9976596617465233]
ic =     [0.355686, 0.354263, 0.352855, 0.351462, 0.350083, 0.348718, 0.347368, 0.346031, 0.344708, 0.343398]



vals = reverse(-Array(1:1:10.0) ./ 10.0 * 2.0 .+ 0.2 / 2.0)

ic_spline = Spline1D(vals, ic)
Y.soil.ϑ_l = ic_spline.(cds.soil.z)
update_aux! = make_update_aux(land)
update_aux!(p, Y, 0.0)

#sim
t0 = FT(0);
N_days = 120
tf = FT(3600 * 24 * N_days)
dt = FT(1);

sv = SavedValues(FT, ClimaCore.Fields.FieldVector)
daily = Array(2:(3600 * 24):(N_days * 3600 * 24))
cb =
    SavingCallback((u, t, integrator) -> copy(integrator.p), sv; saveat = daily)
prob = ODEProblem(ode!, Y, (t0, tf), p);
sol = solve(prob, RK4(), dt = dt, callback = cb); # 8 seconds for 180 days


#plots
ϕ_stem = [
    (Roots.θ_to_p.(sol.u[k].vegetation.θ[1]) .+ 0 * 9800) for
    k in 1:1:length(sol.t)
];
ϕ_leaf = [
    (Roots.θ_to_p.(sol.u[k].vegetation.θ[2]) .+ 18.5 .* 9800) for
    k in 1:1:length(sol.t)
];

lwp_leaf =
    [Roots.θ_to_p.(sol.u[k].vegetation.θ[2]) ./ 1e6 for k in 1:1:length(sol.t)];


    plot1 = plot(seconds ./ 3600 ./ 24, our_year_lwp, label = "LWP (Natan)")
    #plot!(sol.t ./ 3600 ./ 24, ϕ_stem ./ 1e6, label = "ϕ stem, MPa")
    plot!(sol.t ./ 3600 ./ 24, lwp_leaf, label = "LWP (clima) ")
    #plot!(    sol.t ./ 3600 ./ 24,ϕ_leaf ./ 1e6        label = "ϕ leaf, MPa",)
    plot!(
        xlim = [0, maximum(sol.t) ./ 3600 ./ 24],
        ylim = [-3,0],
        xlabel = "t (days since Jan 1)",
        ylabel = "ϕ(MPa)",
    )


    plot!(legend = :bottomright)

    plot2 = plot(
        sol.t ./ 3600 ./ 24,
        [parent(sol.u[k].soil.ϑ_l)[end] for k in 1:1:length(sol.t)],
        label = "10cm",
        xtickfontsize = 5,
        ytickfontsize = 5,
        xlim = [0, maximum(sol.t) ./ 3600 ./ 24],
    )
    plot!(
        sol.t ./ 3600 ./ 24,
        [parent(sol.u[k].soil.ϑ_l)[end - 1] for k in 1:1:length(sol.t)],
        label = "30cm",
    )
#    plot!(
#        sol.t ./ 3600 ./ 24,
#        [parent(sol.u[k].soil.ϑ_l)[end - 2] for k in 1:1:length(sol.t)],
#        label = "50cm",
#    )
 #   plot!(
 #       sol.t ./ 3600 ./ 24,
 #       [parent(sol.u[k].soil.ϑ_l)[end - 3] for k in 1:1:length(sol.t)],
 #       label = "70cm",
 #   )
    plot!(
        sol.t ./ 3600 ./ 24,
        [mean(parent(sol.u[k].soil.ϑ_l)) for k in 1:1:length(sol.t)],
        label = "mean",
    )
    plot!(seconds ./ 3600 ./24, our_year_swc_surface, label = "Natan, surface")
    plot!(seconds ./ 3600 ./24, our_year_swc_obs_surface, label = "obs, surface")
    plot!(seconds ./ 3600 ./24, our_year_swc_column, label = "Natan, mean")
    plot!(legend = :topright, xlabel = "t (days)", ylabel = "θ soil")
    plot!(ylim = [0.0, 0.6])

    plot3 = plot(
        sol.t / 3600 ./ 24,
        precip_function.(sol.t) * 3600 * 39.37,
        xlim = [0, maximum(sol.t) ./ 3600 ./ 24],
        label = "",
    )
    plot!(xlabel = "t (days)", ylabel = "Precip(in/hr)")
    plot4 = plot(
        sol.t / 3600 ./ 24,
        transpiration_function.(sol.t) * 3600 * 39.37,
        xlim = [0, maximum(sol.t) ./ 3600 ./ 24],
        label = "",
    )
    plot!(xlabel = "t (days)", ylabel = "ET(in/hr)")
    plot(plot1, plot2)#, plot3, plot4, layout = 4)
    savefig("./test/LSM/ozark_2005.png")


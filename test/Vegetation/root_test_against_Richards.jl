using DifferentialEquations
using OrdinaryDiffEq:
    ODEProblem,
    solve,
    Rosenbrock23,
    KenCarp4,
    ImplicitEuler
using ModelingToolkit
using BenchmarkTools
using SparsityDetection, SparseArrays
using Plots
#using LaTeXStrings
using Measures
using NLsolve
using Interpolations
using Dierckx

const FT = Float64

function effective_saturation(ν_eff::FT, ϑ_l::FT, θr::FT) where {FT}
    ϑ_l_safe = max(ϑ_l, θr + eps(FT))
    S_l::FT = (ϑ_l_safe - θr) / (ν_eff - θr) # S_l can be > 1
    return S_l # units of [m3 m-3]
end

function volumetric_liquid_fraction(ϑ_l::FT, ν_eff::FT) where {FT} # Do we ever use this, why is it useful in the soil model? I don't think we ever use it. Maybe remove? Afraid of applying pressure correction twice if we use this and the elseif statement for S_l>1 in the absolute_pressure function
    if ϑ_l < ν_eff
        θ_l = ϑ_l
    else
        θ_l = ν_eff
    end
    return θ_l # units of [m3 m-3] 
end

function relative_volumetric_liquid_fraction(θ_l::FT, ν_eff::FT) where {FT}
    rwc = θ_l / ν_eff
    return rwc # units of [m3 m-3]
end

function van_genuchten_volume_to_pressure(α::FT, n::FT, m::FT, S_l::FT) where {FT}
    p = -((S_l^(-FT(1) / m) - FT(1)) * α^(-n))^(FT(1) / n) 
    return p # units of [m]
end

function linear_volume_to_pressure(S_l::FT)      
    p = (S_l - 1) * 5  
    return p # units of [m]
end

function christoffersen_volume_to_pressure_stem(
    θ_l::FT, 
    ν_eff::FT, 
    pi_0::FT, 
    rwc_r::FT, 
    rwc_tlp::FT, 
    rwc_ft::FT, 
    ϵ::FT)  
    rwc = relative_volumetric_liquid_fraction(θ_l, ν_eff) # Christoffersen use rwc ((RWC; g H2 O g−1 H2 O at saturation)). I wonder if we could use S_l instead
    if rwc >= rwc_ft && rwc <= 1    
        p = p_0 - m_cap*(1-rwc)
    elseif rwc >= rwc_tlp && rwc < rwc_ft
        p_sol = -abs(pi_0)*(rwc_ft-rwc_r)/(rwc-rwc_r)
        p_p = abs(pi_0)-epsilon*(rwc_ft-rwc)/(rwc_ft-rwc_r)
        p = p_sol + p_p
    elseif rwc >= rwc_r && rwc < rwc_tlp
        p_sol = -abs(pi_0)*(rwc_ft-rwc_r)/(rwc-rwc_r)
        p = p_sol
    else
        @show("rwc out of bounds")
    end 
    return p # units of [m]
end

function christoffersen_volume_to_pressure_leaf(
    θ_l::FT, 
    ν_eff::FT,
    pi_0::FT, 
    rwc_r::FT, 
    rwc_tlp::FT, 
    ϵ::FT)   
    rwc = relative_volumetric_liquid_fraction(θ_l, ν_eff)
    rwc_ft = FT(1)
    if rwc >= rwc_tlp && rwc <= 1
        p_sol = -abs(pi_0)*(rwc_ft-rwc_r)/(rwc-rwc_r)
        p_p = abs(pi_0)-epsilon*(rwc_ft-rwc)/(rwc_ft-rwc_r)
        p = p_sol + p_p
    elseif rwc >= rwc_r && rwc < rwc_tlp
        p_sol = -abs(pi_0)*(rwc_ft-rwc_r)/(rwc-rwc_r)
        p = p_sol
    else
        @show("rwc out of bounds")
    end 
    return p # units of [m]
end

function absolute_pressure(
    α::FT,
    n::FT,
    m::FT,
    θ_r::FT,
    ϑ_l::FT,
    ν_eff::FT,
    S_s::FT
) where {FT}
    S_l = effective_saturation(ν_eff, ϑ_l, θ_r)
    if S_l <= FT(1.0)
        p = van_genuchten_volume_to_pressure(α, n, m, S_l)
    else
        p = (ϑ_l - ν_eff) / S_s
    end
    return p # units of [m]
end

function absolute_pressure(
    ϑ_l::FT,
    θ_r::FT,
    ν_eff::FT,
    S_s::FT
) where {FT}
    S_l = effective_saturation(ν_eff, ϑ_l, θ_r)
    if S_l <= FT(1.0)
        p = linear_volume_to_pressure(S_l)
    else
        p = (ϑ_l - ν_eff) / S_s
    end
    return p # units of [m]
end

#=
function absolute_pressure(
    θ_l::FT, 
    ϑ_l::FT,
    ν_eff::FT, 
    pi_0::FT, 
    rwc_r::FT, 
    rwc_tlp::FT, 
    rwc_ft::FT, 
    ϵ::FT
) where {FT}
    S_l = effective_saturation(ν_eff, ϑ_l, θ_r)
    if S_l <= FT(1.0)
        p = christoffersen_volume_to_pressure_stem(θ_l, ν_eff, pi_0, rwc_r, rwc_tlp, rwc_ft, ϵ)
    else
        p = (ϑ_l - ν_eff) / S_s
    end
    return p # units of [m]
end

function absolute_pressure(
    θ_l::FT, 
    ϑ_l::FT,
    ν_eff::FT, 
    pi_0::FT, 
    rwc_r::FT, 
    rwc_tlp::FT, 
    rwc_ft::FT, 
    ϵ::FT
) where {FT}
    S_l = effective_saturation(ν_eff, ϑ_l, θ_r)
    if S_l <= FT(1.0)
        p = christoffersen_volume_to_pressure_leaf(θ_l, ν_eff, pi_0, rwc_r, rwc_tlp, ϵ)
    else
        p = (ϑ_l - ν_eff) / S_s
    end
    return p # units of [m]
end
=#

function van_genuchten_hydraulic_conductivity(
    Ksat::FT,
    m::FT,
    S_l::FT
) where {FT}
    if S_l < FT(1)
        K = sqrt(S_l) * (FT(1) - (FT(1) - S_l^(FT(1) / m))^m)^FT(2)
    else
        K = FT(1)
    end
    return K * Ksat # units of [m s-1]
end

function modified_logistic_hydraulic_conductivity(
    Ksat::FT,
    a::FT,
    b::FT,
    S_l::FT
) where {FT}
    if S_l < FT(1)
        K = (a + 1)/a * (FT(1) - FT(1)/(FT(1) + a * exp(b * p)))
    else
        K = FT(1)
    end
    return K * Ksat # units of [m s-1]
end

function rhs_flux!(dY, Y, paramset1, t)
    ν,
    vg_α,
    vg_n,
    vg_m,
    Ksat,
    S_s,
    θ_r,
    top_flux_bc,
    bot_flux_bc = paramset1
    S_l = effective_saturation.(ν, Y, θ_r)
    K = van_genuchten_hydraulic_conductivity.(Ksat,vg_m,S_l) # units of [m s-1]
    p = absolute_pressure.(vg_α,vg_n,vg_m,θ_r,Y,ν,S_s) .+ z  # units [m]

    @inbounds for i in 1:1:n
        ip1, im1 = i+1, i-1
        Fi_ph::FT = ip1 == n+1 ? FT.(top_flux_bc) : FT(-2.0)/(FT(1.0)/K[ip1]+FT(1.0)/K[i])*(p[ip1]-p[i])/Δz
        Fi_mh::FT = im1  == 0 ? FT.(bot_flux_bc) : FT(-2.0)/(FT(1.0)/K[i]+FT(1.0)/K[im1])*(p[i]-p[im1])/Δz
        dY[i] = -FT(1.0)/Δz*(Fi_ph-Fi_mh)
    end
end

function flux(
    z1::FT,
    z2::FT, 
    p1::FT, 
    p2::FT, 
    a1::FT, 
    a2::FT, 
    b1::FT, 
    b2::FT, 
    K_sat1::FT, 
    K_sat2::FT) where {FT}  
        u1 = a1 * exp(b1 * p1) 
        u2 = a1 * exp(b1 * p2) 
        num1 = log(u1 + FT(1))
        num2 = log(u2 + FT(1))
        c1 = K_sat1 * (a1 + FT(1)) / a1
        term1 = -c1 / b1 * (num2 - num1) /(z2 - z1)

        c2 = K_sat2 * (a2 + FT(1)) / a2
        term2_up = -c2 * (FT(1) - FT(1) / (FT(1) + a2*exp(b2 * p2)))
        term2_do = -c1 * (FT(1) - FT(1) / (FT(1) + a1*exp(b1 * p1)))
        term2 = (term2_up + term2_do)/2
        flux = term1 + term2  
    return flux  # units of [m s-1]
end

function flux_Yujie(
    z1::FT,
    z2::FT,
    p1::FT,
    p2::FT,
    a1::FT,
    b1::FT,
    K_sat1::FT
)::FT where {FT}
    u1, u2, A, B, flux_approx = vc_integral_approx(z1, z2, p1, p2, a, b, K_sat1)
    @show(u1)
    @show(u2)
    @show(A)
    @show(B)
    @show(flux_approx)
    flux = vc_integral(u1, u2, A, B, flux_approx)
    return flux  # units of [m s-1]
end

function vc_integral_approx(
    z1::FT,
    z2::FT,
    p1::FT,
    p2::FT,
    a::FT,
    b::FT,
    K_sat::FT,
) where {FT}
    u1 = a * exp(b * p1)
    u2 = a * exp(b * p2)
    num1 = log(u1 + FT(1))
    num2 = log(u2 + FT(1))
    c = K_sat * (a + FT(1)) / a
    d = z2 - z1
    flux_approx = -c / (b * d) * (num2 - num1) * (p2 - p1 + d) / (p2 - p1) # this is NaN if p2 = p1
    A = c + flux_approx
    B = -c * flux_approx / (b * d * A)
    @show(u1)
    @show(u2)
    @show(A)
    @show(B)
    @show(flux_approx)
    return u1, u2, A, B, flux_approx  # units of [m s-1]
end

function vc_integral(u1::FT, u2::FT, A::FT, B::FT, flux_approx::FT) where {FT}
    flux = B * log((u2 * A + flux_approx) / (u1 * A + flux_approx))
    return flux # units of [m-1]
end

function roots(dy,y,paramset2,t)
    p = absolute_pressure.(vg_α,vg_n,vg_m,θ_r,y,ν,S_s)
    m = Int64(paramset2[end-2])
    z2 = paramset2[end-1]
    z3 = paramset2[end]

    flux_stem = ones(m+1,1)
    flux_stem[1] = bot_flux_bc
    flux_stem[end] = top_flux_bc
    
    for i in 1:m-1
        flux_stem[i+1] = flux(z2[i], z2[i+1], p[i], p[i+1], a, a, b, b, Ksat, Ksat)
    end
    
    # Assuming uniform conducting areas so not including area term
    for j in 1:m
        dy[j] = 1/(z3[j+1] - z3[j]) * (flux_stem[j] - flux_stem[j+1])
    end  
end

function roots2(dy,y,paramset2,t)
    p = absolute_pressure.(vg_α,vg_n,vg_m,θ_r,y,ν,S_s)
    m = Int64(paramset2[end-1])
    z2 = paramset2[end]

    flux_stem = ones(m+1,1)
    flux_stem[1] = bot_flux_bc
    flux_stem[end] = top_flux_bc
    
    for i in 1:m-1
        flux_stem[i+1] = flux(z2[i], z2[i+1], p[i], p[i+1], a, a, b, b, Ksat, Ksat)
    end
    
    # Assuming uniform conducting areas so not including area term
    for j in 1:m
        dy[j] = 1/(z2[j+1] - z2[j]) * (flux_stem[j] - flux_stem[j+1])
    end  
end

function roots_Yujie(dy,y,paramset2,t)
    p = absolute_pressure.(vg_α,vg_n,vg_m,θ_r,y,ν,S_s)
    m = Int64(paramset2[end-1])
    z2 = paramset2[end]

    flux_stem = ones(m+2,1)
    flux_stem[1] = bot_flux_bc
    flux_stem[end] = top_flux_bc
    
    for i in 1:m-1
        flux_stem[i+1] = flux_Yujie(z2[i], z2[i+1], p[i], p[i+1], a, b, Ksat)
    end

    # Assuming uniform conducting areas so not including area term
    for j in 1:m
        dy[j] = 1/(z2[j+1] - z2[j]) * (flux_stem[j] - flux_stem[j+1])
    end   
end

# Model parameters
ν = FT(0.495)
Ksat = FT(0.0443 / 3600 / 100) # units of [m s-1]
S_s = FT(1e-3) # units of [m-1]
vg_n = FT(2.0)
vg_α = FT(2.6) # units of [m-1]
vg_m = FT(1) - FT(1) / vg_n
a = FT(0.3968) # modified logistic function parameters fitted to VG, for plant hydraulics model
b = FT(7.692)
θ_r = FT(0)
zmax = FT(10)
zmin = FT(0)

const day = (60^2) * 24
const hour = 60^2
const n = 60
const Δz = (zmax-zmin)/n
const z = Array((zmin+Δz/FT(2.0)):Δz:(zmax-Δz/FT(2.0)))

top_flux_bc = FT(0.0)
bot_flux_bc = FT(0.0)

paramset1 = [ν, vg_α, vg_n, vg_m, Ksat, S_s, θ_r, top_flux_bc, bot_flux_bc]

t0 = FT(0)
nday = FT(1)
tf = FT(60 * 60 * 24 * nday)
dt = FT(10)
ν_vec = FT(0.495) .+ zeros(n)
θ0 =  FT(0.494) .+ zeros(n)

prob = ODEProblem(rhs_flux!, θ0, (t0, tf), paramset1);
i = Array(1:1:n)
iu = Array(1:1:n-1)
il = Array(2:1:n)
i1 = vcat(i,iu,il)
i2 = vcat(i,il,iu)
jac_sparsity = sparse(i1,i2, 1.0)
f = ODEFunction(rhs_flux!;jac_prototype=jac_sparsity);
prob_sparse = ODEProblem(f, θ0, (t0,tf), paramset1);

# Solve
##### Truth against model for varying number of compartments, and fixed simulation time
@time begin
    truth = solve(prob,dt = dt,SSPRK33(),save_every_step=false)
    end
    
RE_against_PHM = plot(z,truth.u[1],label="t0",xlabel="Stem height [m]",ylabel="ϑ stem [m3 m-3]",legend=:bottomleft,dpi=300)
plot!(RE_against_PHM,z,truth.u[end],label=string("truth n=", string(n)))

# Simulation length
#tspan = (0.0,tf)
alg = Euler()
alg_name = "Euler"

num_compartments = Array([2,10,20,40,60])
absolute_error_as_function_of_height=plot(xlabel="Stem height [m]",ylabel="ϑ truth - ϑ model [m3 m-3]")

for i in 1:length(num_compartments)
m = num_compartments[i] # number of compartments
y0 =  FT(0.494) .+ zeros(m)
Δz2 = (zmax-zmin)/m # size of compartments

# To use with "Flux2" function
z2 = Array(zmin:Δz2:zmax)
z3 = Array((zmin+Δz2/FT(2.0)):Δz2:(zmax-Δz2/FT(2.0)))

# To use with "Flux" function
#z2 = Array((zmin+Δz2/FT(2.0)):Δz2:(zmax-Δz2/FT(2.0))) # position of pressure values
#z3 = Array(zmin:Δz2:zmax) # position of boundaries

# Solve the problem
paramset2 = [ν, vg_α, vg_n, vg_m, Ksat, S_s, θ_r, top_flux_bc, bot_flux_bc, a, b, m, z2] #, z3] Use z3 with "Flux" function
@time begin
sol = solve(ODEProblem(roots2,y0,tf,paramset2),alg,adaptive=false,dt=dt)
end

theta = ones(m,length(sol.u))
for j in 1:m
    theta[j,:] = reduce(hcat,sol.u)[j,:]
end

theta_start=theta[:,1]
theta_end=theta[:,end] 

# Spline parameters
if m>7
    k=5
else
    k=m-1
end

spl = Spline1D(z3, theta_end, k=k)
theta_interpolated(zs) = spl(zs)
A = [theta_interpolated(zs) for zs in z3]
interp_linear = LinearInterpolation(z3, A)
extrap = LinearInterpolation(z3, A,extrapolation_bc=Line()) 
extrap_theta_at_zmin= extrap(z[1]) 
extrap_theta_at_zmax= extrap(z[end])
z4 = [z[1]; z3; z[end]]
theta_end = [extrap_theta_at_zmin; theta_end; extrap_theta_at_zmax]

plot!(RE_against_PHM,z4,theta_end,line=(:dash,1),label=string("model n=", string(m)),dpi=300)
scatter!([z4[1],z4[2],z4[end-1],z4[end]],[theta_end[1],theta_end[2],theta_end[end-1],theta_end[end]],label="",markershape = :circle, markersize = 1.5, markercolor = :black)

# Error as a function of height
theta_end_at_z = extrap(z)
abs_error = truth.u[end] - theta_end_at_z
plot!(absolute_error_as_function_of_height,z,abs_error,label=string("n=",string(m)),legend=:bottomleft,dpi=300)
end

plot!(RE_against_PHM,z,ν_vec,label="porosity",xlabel="Stem height [m]", ylabel="ϑ [m3 m-3]")
savefig(RE_against_PHM,string("RE_against_PHM_for_varying_num_compartments_", string(n), "_elements_sim_time_", string(nday), "_days.png"))
savefig(absolute_error_as_function_of_height,string("absolute_error_as_function_of_height_", string(n), "_elements_sim_time_", string(nday), "_days.png"))

################### 
#Model against truth for fixed number of compartments, and varying simulation time

m = 60
Δz2 = (zmax-zmin)/m # size of compartments
z2 = Array((zmin+Δz2/FT(2.0)):Δz2:(zmax-Δz2/FT(2.0))) # position of pressure values
z3 = Array(zmin:Δz2:zmax) # position of boundaries
y0 =  FT(0.494) .+ zeros(m)

# Solve the problem
paramset2 = [ν, vg_α, vg_n, vg_m, Ksat, S_s, θ_r, top_flux_bc, bot_flux_bc, a, b, m, z2, z3]
@time begin
sol = solve(ODEProblem(roots,y0,tf,paramset2),alg,adaptive=false,dt=dt)
end

theta = ones(m,length(sol.u))
for j in 1:m
    theta[j,:] = reduce(hcat,sol.u)[j,:]
end

###### Model vs truth as function of time, over 6 hours
label_truth=["truth 0h";"truth 1h";"truth 2h"; "truth 3h"; "truth 4h"; "truth 5h"; "truth 6h"] 
label_model=["model 0h";"model 1h";"model 2h"; "model 3h"; "model 4h"; "model 5h"; "model 6h"] 
RE_against_PHM_6h = plot()
let k=0
    for i in Int64(1):Int64(hour/dt):Int64(7*hour/dt) 
        k=k+1
        plot!(RE_against_PHM_6h,z,truth.u[i],label=label_truth[k], legendfontsize=10, color=RGB((k+2)/10, 0.3, 0.4))
        plot!(RE_against_PHM_6h,z2,theta[:,i],line=(:dash, 1),label=label_model[k],color=RGB((k+2)/10, 0.3, 0.4),legend=:bottomleft,dpi=300,xlabel="Stem height [m]",ylabel="Augmented liquid fraction in stem [m3 m-3]") #line=(:dot, 1)
    end
end
plot!(RE_against_PHM_6h,z,ν_vec,label="porosity")
savefig(RE_against_PHM_6h,string("RE_against_PHM_6h_",string(m),"_compartments_",string(n),"elements.png"))

###### Model vs truth as function of time, over 20 days
#=
label_truth=["truth 1 day";"truth 5 days";"truth 10 days";"truth 15 days";"truth 20 days"] 
label_model=["model 1 day";"model 5 days";"model 10 days";"model 15 days";"model 20 days"] 
RE_against_PHM_20_days = plot()
let k=0
    for i in Int64(day/dt):Int64(5*day/dt):Int64(25*day/dt) 
        k=k+1
        plot!(RE_against_PHM_20_days,z,truth.u[i],label=label_truth[k], legendfontsize=10, color=RGB((k+2)/10, 0.3, 0.4))
        plot!(RE_against_PHM_20_days,z2,theta[:,i],line=(:dash, 1),label=label_model[k],color=RGB((k+2)/10, 0.3, 0.4),legend=:bottomleft,dpi=300,xlabel="Stem height [m]",ylabel="Augmented liquid fraction in stem [m3 m-3]")
    end
end
plot!(RE_against_PHM_20_days,z,ν_vec,label="porosity")
savefig(RE_against_PHM_20_days,"RE_against_PHM_20_days_60_compartments_100_elements.png")
=#

# Useful functions
# to annotate figures: annotate!(RE_against_PHM_6h, 0.1, 0.4910, text("full line: truth", :left, 10))
# to concatenate: string()
# to us Latex font: absolute_error_as_function_of_height=plot(xlabel=L"\mathrm{ Stem height [m] }",ylabel=L"\mathrm{ \theta_{truth} - \theta_{model} [m^3 m^{-3}] }")
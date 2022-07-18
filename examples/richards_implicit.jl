using LinearAlgebra
using OrdinaryDiffEq

using ClimaCore: Operators, Geometry, Spaces, Fields
using ClimaCore.Utilities: half
if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
using ClimaLSM
using ClimaLSM.Soil
using ClimaLSM.Domains: HybridBox, Column
include("ordinary_diff_eq_bug_fixes.jl")
include("autodiff_test.jl")
using UnPack
const compose = Operators.ComposeStencils()
const apply = Operators.ApplyStencil()
# Allow one() to be called on vectors.
Base.one(::T) where {T <: Geometry.AxisTensor} = one(T)
Base.one(::Type{T}) where {T′, A, S, T <: Geometry.AxisTensor{T′, 1, A, S}} =
    T(axes(T), S(one(T′)))

# As long as the boundary conditions are independent of state, we
# don't care about their values in computing Operator2Stencil
FT = Float64
gradc2f_op = Operators.GradientC2F(
    bottom = Operators.SetValue(FT(0.0)),
    top = Operators.SetValue(FT(0.0)),
)
gradc2f_stencil = Operators.Operator2Stencil(gradc2f_op)
divf2c_op = Operators.DivergenceF2C(
    top = Operators.SetValue(Geometry.Covariant3Vector(FT(0))),
    bottom = Operators.SetValue(Geometry.Covariant3Vector(FT(0))),
)
divf2c_stencil = Operators.Operator2Stencil(divf2c_op)
interpc2f_op = Operators.InterpolateC2F(
    bottom = Operators.Extrapolate(),
    top = Operators.Extrapolate(),
)
interpc2f_stencil = Operators.Operator2Stencil(interpc2f_op)

struct TridiagonalJacobian{R, J1, J2}
    # reference to dtγ, which is specified by the ODE solver
    dtγ_ref::R
    ∂ϑt∂ϑ::J1
    ∂ϑt∂ϑ_column_array::J2
    flag::Symbol
    test_flag::Bool
end

function TridiagonalJacobian(Y, flag::Symbol; test_flag::Bool = false)
    FT = eltype(Y.soil.ϑ_l)
    dtγ_ref = Ref(zero(FT))
    space = axes(Y.soil.ϑ_l)
    N = Spaces.nlevels(space)
    Jtype = Operators.StencilCoefs{-1, 1, NTuple{3, FT}}
    J = Fields.Field(Jtype, space)
    column_array = Tridiagonal(
        Array{FT}(undef, N - 1),
        Array{FT}(undef, N),
        Array{FT}(undef, N - 1),
    )
    args = (dtγ_ref, J, column_array)
    return TridiagonalJacobian{typeof.(args)...}(args..., flag, test_flag)
end
Base.similar(w::TridiagonalJacobian) = w
to_scalar_coefs(vector_coefs) = map(vector_coef -> vector_coef.u₃, vector_coefs)
function make_Wfact(model::RichardsModel)
    function Wfact!(W, Y, p, dtγ, t)
        (; dtγ_ref, ∂ϑt∂ϑ, flag, test_flag) = W
        dtγ_ref[] = dtγ
        # ∂ϑ∂t = ϑt = div(I(K(θ)) grad(h(θ)))
        if flag == :exact
            @. ∂ϑt∂ϑ = compose(
                divf2c_stencil(
                    one.(gradc2f_op.(p.soil.ψ + model.coordinates.z)),
                ),
                to_scalar_coefs(
                    interpc2f_op(dKdθ(Y.soil.ϑ_l, ν, θ_r, vg_m, K_sat)) *
                    interpc2f_stencil(one(Y.soil.ϑ_l)) *
                    gradc2f_op(p.soil.ψ + model.coordinates.z) +
                    interpc2f_op(
                        p.soil.K *
                        dψdθ(Y.soil.ϑ_l, ν, θ_r, vg_α, vg_n, vg_m, S_s),
                    ) * gradc2f_stencil(one.(Y.soil.ϑ_l)),
                ),
            )
        elseif flag == :modifiedpicard
            @. ∂ϑt∂ϑ = compose(
                divf2c_stencil(
                    one.(gradc2f_op.(p.soil.ψ + model.coordinates.z)),
                ),
                to_scalar_coefs(
                    interpc2f_op(
                        p.soil.K *
                        dψdθ(Y.soil.ϑ_l, ν, θ_r, vg_α, vg_n, vg_m, S_s),
                    ) * gradc2f_stencil(one.(Y.soil.ϑ_l)),
                ),
            )
        else
            println("Unsupported Jacobian approximation")
        end
        i, j, h = 1, 1, 1

        if test_flag
            args = (implicit_tendency!, Y, p, t, i, j, h)
            display(matrix_column(∂ϑt∂ϑ, axes(Y.soil.ϑ_l), i, j, h))
            display(
                exact_column_jacobian_block(
                    args...,
                    (:soil, :ϑ_l),
                    (:soil, :ϑ_l),
                ),
            )

            @assert matrix_column(∂ϑt∂ϑ, axes(Y.soil.ϑ_l), i, j, h) ≈
                    exact_column_jacobian_block(
                args...,
                (:soil, :ϑ_l),
                (:soil, :ϑ_l),
            )
        end

        # CLM approach
        # Newton's method with one iteration, and use the previous time step as your first guess.

    end
    return Wfact!
end

function dKdθ(θ, ν, θ_r, vg_m, Ksat)
    S = (θ - θ_r) / (ν - θ_r)
    if S < 1
        f = 1.0 - (1.0 - S^(1.0 / vg_m))^vg_m
        f1 = f^2.0 / 2.0 / S^0.5
        f2 = 2 * S^(1 / vg_m - 1 / 2) * f / (1 - S^(1 / vg_m))^(1.0 - vg_m)
        return (f1 + f2) * Ksat / (ν - θ_r)
    else
        return zero(typeof(S))
    end

end
#    if S < 1
#        return (((1 - (1 - S^(1 / vg_m))^vg_m))^2 / 2 / S^0.5
#                + 2 * S^(1 / vg_m - 1 / 2) * (1 - (1 - S^(1 / vg_m))^vg_m) / (1 - S^(1 / vg_m))^(1 - vg_m)
#                )* Ksat / (ν - θ_r)
#    else
#        return zero(typeof(S))
#    end

#end

function dψdθ(θ, ν, θ_r, vg_α, vg_n, vg_m, S_s)
    S = (θ - θ_r) / (ν - θ_r)
    if S < 1.0
        return 1.0 / (vg_α * vg_m * vg_n) / (ν - θ_r) *
               (S^(-1 / vg_m) - 1)^(1 / vg_n - 1) *
               S^(-1 / vg_m - 1)
    else
        return 1.0 / S_s
    end
end

#    if S < 1
#        return 1 / (vg_α * vg_m * vg_n) / (ν - θ_r) *
#            (S^(-1 / vg_m) - 1)^(1 / vg_n - 1) *
#            S^(-1 / vg_m - 1)
#    else
#        return 1 / S_s
#    end
#end


function make_implicit_tendency(model::RichardsModel)
    function implicit_tendency!(dY, Y, p, t)
        @unpack ν, vg_α, vg_n, vg_m, K_sat, S_s, θ_r = model.parameters
        # when does aux get updated?
        (; K, ψ) = p.soil
        if eltype(Y) <: Dual
            K = similar(Y.soil.ϑ_l)
            ψ = similar(Y.soil.ϑ_l)
        end
        @. K = hydraulic_conductivity(
            K_sat,
            vg_m,
            effective_saturation(ν, Y.soil.ϑ_l, θ_r),
        )
        @. ψ = pressure_head(vg_α, vg_n, vg_m, θ_r, Y.soil.ϑ_l, ν, S_s)

        top_flux_bc, bot_flux_bc =
            boundary_fluxes(model.boundary_conditions, p, t)
        z = model.coordinates.z
        interpc2f = Operators.InterpolateC2F()
        gradc2f_water = Operators.GradientC2F()
        divf2c_water = Operators.DivergenceF2C(
            top = Operators.SetValue(Geometry.WVector(top_flux_bc)),
            bottom = Operators.SetValue(Geometry.WVector(bot_flux_bc)),
        )
        @. dY.soil.ϑ_l = -(divf2c_water(-interpc2f(K) * gradc2f_water(ψ + z)))
    end
    return implicit_tendency!
end


function make_explicit_tendency(model::RichardsModel)
    function explicit_tendency!(dY, Y, p, t)
        #when does aux get updated?
        @. p.soil.K = hydraulic_conductivity(
            K_sat,
            vg_m,
            effective_saturation(ν, Y.soil.ϑ_l, θ_r),
        )
        @. p.soil.ψ = pressure_head(vg_α, vg_n, vg_m, θ_r, Y.soil.ϑ_l, ν, S_s)
        hdiv = Operators.WeakDivergence()
        hgrad = Operators.Gradient()
        @. dY.soil.ϑ_l +=
            -hdiv(-p.soil.K * hgrad(p.soil.ψ + model.coordinates.z))
        Spaces.weighted_dss!(dY.soil.ϑ_l)
    end
    return explicit_tendency!
end

function linsolve!(::Type{Val{:init}}, f, u0; kwargs...)
    function _linsolve!(x, A, b, update_matrix = false; kwargs...)
        # TODO: Do this with stencil_solve!.
        ϑ_l = x.soil.ϑ_l
        (; ∂ϑt∂ϑ, ∂ϑt∂ϑ_column_array) = A
        Ni, Nj, _, _, Nh = size(Spaces.local_geometry_data(axes(ϑ_l)))
        for h in 1:Nh, j in 1:Nj, i in 1:Ni
            ϑ_column_view = parent(Spaces.column(ϑ_l, i, j, h))
            ∂ϑt∂ϑ_column = Spaces.column(∂ϑt∂ϑ, i, j, h)
            @views ∂ϑt∂ϑ_column_array.dl .= parent(∂ϑt∂ϑ_column.coefs.:1)[2:end]
            ∂ϑt∂ϑ_column_array.d .= parent(∂ϑt∂ϑ_column.coefs.:2)
            @views ∂ϑt∂ϑ_column_array.du .=
                parent(∂ϑt∂ϑ_column.coefs.:3)[1:(end - 1)]
            ldiv!(lu!(∂ϑt∂ϑ_column_array), ϑ_column_view)
        end
    end
end



domain = HybridBox(;
    zlim = (-10.0, 0.0),
    xlim = (0.0, 100.0),
    ylim = (0.0, 100.0),
    nelements = (10, 10, 10),
    npolynomial = 1,
    periodic = (true, true),
)
ν = FT(0.495)
K_sat = FT(0.0443 / 3600 / 100) # m/s
S_s = FT(1e-3) #inverse meters
vg_n = FT(2.0)
vg_α = FT(2.6) # inverse meters
vg_m = FT(1) - FT(1) / vg_n
θ_r = FT(0)

soil_domain = domain
top_flux_bc = FT(0.0)
bot_flux_bc = FT(0.0)
sources = ()
boundary_conditions = FluxBC{FT}(top_flux_bc, bot_flux_bc)
params = RichardsParameters{FT}(ν, vg_α, vg_n, vg_m, K_sat, S_s, θ_r)

soil = RichardsModel{FT}(;
    parameters = params,
    domain = soil_domain,
    boundary_conditions = boundary_conditions,
    sources = sources,
)

Y, p, coords = initialize(soil)
function init_soil!(Y, z, params)
    function hydrostatic_profile(
        z::FT,
        params::RichardsParameters{FT},
    ) where {FT}
        @unpack ν, vg_α, vg_n, vg_m, θ_r = params
        #unsaturated zone only, assumes water table starts at z_∇
        z_∇ = FT(-10)# matches zmin
        S = FT((FT(1) + (vg_α * (z - z_∇))^vg_n)^(-vg_m))
        ϑ_l = S * (ν - θ_r) + θ_r
        return FT(ϑ_l)
    end
    Y.soil.ϑ_l .= hydrostatic_profile.(z, Ref(params))
end

init_soil!(Y, coords.z, soil.parameters)
update_aux! = make_update_aux(soil)
update_aux!(p, Y, 0.0)

W = TridiagonalJacobian(Y, :exact; test_flag = true)
Wfact! = make_Wfact(soil)
t_start = 0.0
t_end = 100.0
dt = 10.0
implicit_tendency! = make_implicit_tendency(soil)
explicit_tendency! = make_explicit_tendency(soil)
jac_kwargs = (; jac_prototype = W, Wfact = Wfact!)
alg_kwargs = (; linsolve = linsolve!)
problem = SplitODEProblem(
    ODEFunction(
        implicit_tendency!;
        jac_kwargs...,
        tgrad = (∂Y∂t, Y, p, t) -> (∂Y∂t .= FT(0)),
    ),
    explicit_tendency!,
    Y,
    (t_start, t_end),
    p,
)
integrator = init(
    problem,
    Rosenbrock23(; alg_kwargs...);
    dt = dt,
    adaptive = false,
    progress = true,
    progress_steps = 1,
)

sol = @timev OrdinaryDiffEq.solve!(integrator)

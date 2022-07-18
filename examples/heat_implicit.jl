using LinearAlgebra

using ClimaCore: Operators, Geometry, Spaces, Fields
using ClimaCore.Utilities: half
using ClimaLSM
using ClimaLSM.Domains: HybridBox, Column
using ClimaCore.Utilities: half

const compose = Operators.ComposeStencils()
const apply = Operators.ApplyStencil()
# Allow one() to be called on vectors.
Base.one(::T) where {T <: Geometry.AxisTensor} = one(T)
Base.one(::Type{T}) where {T′, A, S, T <: Geometry.AxisTensor{T′, 1, A, S}} =
    T(axes(T), S(one(T′)))


function implicit_tendency!(dY,Y,p,t)
    bcs_bottom = Operators.SetValue(Geometry.WVector(FT(0.0)))
    bcs_top = Operators.SetValue(Geometry.WVector(FT(1.0)))
    
    gradc2f = Operators.GradientC2F()
    divf2c = Operators.DivergenceF2C(bottom = bcs_bottom, top = bcs_top)

    return @. dY = divf2c(gradc2f(Y.T))
end


function explicit_tendency!(dY,Y,p,t)
    hdiv = Operators.WeakDivergence()
    hgrad = Operators.Gradient()
    @. dY.T += - hdiv(-hgrad(Y.T))
    Spaces.weighted_dss!(dY.T)
end

# As long as the boundary conditions are independent of state, we
# don't care about their values in computing Operator2Stencil
FT = Float64
grad_op = Operators.GradientC2F(
    bottom = Operators.SetValue(FT(0.0)),
    top = Operators.SetValue(FT(0.0))
)
grad_stencil = Operators.Operator2Stencil(grad_op)
div_op = Operators.DivergenceF2C(
    top = Operators.SetValue(Geometry.Covariant3Vector(FT(0))),
    bottom = Operators.SetValue(Geometry.Covariant3Vector(FT(0))),
)
div_stencil = Operators.Operator2Stencil(div_op)
interpc2f = Operators.InterpolateC2F(
    bottom = Operators.Extrapolate(),
    top = Operators.Extrapolate(),
)

struct TridiagonalJacobian{R, J1}
    # reference to dtγ, which is specified by the ODE solver
    dtγ_ref::R
    ∂Tt∂T::J1
end

function TridiagonalJacobian(Y)
    FT = eltype(Y.T)
    dtγ_ref = Ref(zero(FT))
    space = axes(Y.T)
    Jtype = Operators.StencilCoefs{-1, 1, NTuple{3, FT}}
    J = Fields.Field(Jtype, space)
    args = (dtγ_ref,J)
    return TridiagonalJacobian{typeof.(args)...}(args...)
end
Base.similar(w::TridiagonalJacobian) = w
to_scalar_coefs(vector_coefs) =
    map(vector_coef-> vector_coef.u₃, vector_coefs)
function Wfact!(W, Y, p, dtγ, t)
    (; dtγ_ref, ∂Tt∂T) = W
    dtγ_ref[] = dtγ
    # ∂T∂t = Tt = div(grad(T))
    @. ∂Tt∂t = compose(div_stencil(one.(grad_op.(Y.T))), to_scalar_coefs(grad_stencil(one.(Y.T))))

end

#domain = HybridBox(;
#                   zlim = (0.0,1.0),
#                   xlim = (0.0,100.0),
#                   ylim = (0.0,100.0),
#                   nelements = (10,10,10),
#                   npolynomial = 1,
#                   periodic = (true,true)
#                   )
domain = Column(; zlim = (0.0,1.0), nelements = 10)
cc = ClimaLSM.Domains.coordinates(domain)
T0 = cc.z.^2.0
                   
tspan = (0, 10)
Y = Fields.FieldVector(;:T => T0,);
p = nothing
W = TridiagonalJacobian(Y)
### Need a linsolve...this doesnt do anything yet.
function linsolve!(::Type{Val{:init}}, f, u0; kwargs...)
    function _linsolve!(x, A, b, update_matrix = false; kwargs...)
        (; dtγ_ref, ∂Tt∂T) = W
        A = Tridiagonal(∂Tt∂t.upper_diagonal,
                        ∂Tt∂t.diagonal,
                        ∂Tt∂t.lower_diagonal)
        dtγ = dtγ_ref[]
    end
end
# copied code for running, eventually
#=
jac_kwargs = (; jac_prototype = W, Wfact = Wfact!)
max_newton_iters = 2
alg_kwargs = (; linsolve = linsolve!)
alg_kwargs =
            (; alg_kwargs..., nlsolve = NLNewton(; max_iter = max_newton_iters))
problem = SplitODEProblem(
    ODEFunction(
        implicit_tendency!;
        jac_kwargs...,
        tgrad = (∂Y∂t, Y, p, t) -> (∂Y∂t .= FT(0)),
    ),
    explicit_tendency!,
    Y0,
    tspan,
    p,
)
=#

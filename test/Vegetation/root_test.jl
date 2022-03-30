using Test
using Statistics
using NLsolve
using OrdinaryDiffEq: ODEProblem, solve, RK4
using ClimaCore

if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
using ClimaLSM
using ClimaLSM.Domains: RootDomain
using ClimaLSM.Roots

FT = Float64

@testset "Root model integration tests" begin
    a_root = FT(13192)
    a_stem = FT(515.5605)
    b_root = FT(2.1079 / 1e6) # Inverse Pa
    b_stem = FT(0.9631 / 1e6) # Inverse Pa
    h_leaf = FT(0.01) #10mm, guess
    h_stem = FT(0.5)
    K_max_stem = FT(3.75e-9)
    K_max_root = FT(1.48e-8)
    SAI = FT(0.1)
    RAI = FT(0.1)
    LAI = FT(0.3)
    z_leaf = FT(0.5) # height of leaf
    z_root_depths = FT.([-1.0]) # m, rooting depth
    z_bottom_stem = FT(0.0)

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
        (z) -> 1.0,
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

    p_soil0 = FT.([-2e6])
    transpiration =
        PrescribedTranspiration{FT}((t::FT) -> leaf_transpiration(t))
    root_extraction = PrescribedSoilPressure{FT}((t::FT) -> p_soil0)
    roots = Roots.RootsModel{FT}(;
        domain = root_domain,
        param_set = param_set,
        root_extraction = root_extraction,
        transpiration = transpiration,
    )

    # Set system to equilibrium state by setting LHS of both ODEs to 0


    function f!(F, Y)
        T0 = 1e-5
        flow_in_stem = sum(
            ground_area_flux.(
                z_root_depths,
                z_bottom_stem,
                p_soil0,
                Y[1],
                a_root,
                b_root,
                K_max_root,
                RAI,
            ),
        )
        flow_out_stem = ground_area_flux(
            z_bottom_stem,
            z_leaf,
            Y[1],
            Y[2],
            a_stem,
            b_stem,
            K_max_stem,
            SAI,
        )
        F[1] = flow_in_stem - flow_out_stem
        F[2] = flow_out_stem - T0 * LAI
    end

    soln = nlsolve(f!, [-1e6, -0.9e6], ftol = 1e-10)
    p_stem_ini = soln.zero[1]
    p_leaf_ini = soln.zero[2]

    theta_stem_ini = p_to_θ(p_stem_ini)
    theta_leaf_ini = p_to_θ(p_leaf_ini)
    y0 = FT.([theta_stem_ini, theta_leaf_ini])
    Y, p, coords = initialize(roots)
    Y.vegetation.θ .= y0

    root_ode! = make_ode_function(roots)

    t0 = FT(0)
    tf = FT(1200)
    dt = FT(1)

    prob = ODEProblem(root_ode!, Y, (t0, tf), p)
    sol = solve(prob, RK4(), dt = dt)

    dY = similar(Y)
    root_ode!(dY, Y, p, 0.0)
    @test sqrt(mean(dY.vegetation.θ .^ 2.0)) < 1e-8 # starts in equilibrium


    y_1 = reduce(hcat, sol.u)[1, :]
    y_2 = reduce(hcat, sol.u)[2, :]
    p_stem = θ_to_p.(y_1)
    p_leaf = θ_to_p.(y_2)

    function f2!(F, Y)
        p_soilf = p_soil0
        Tf = 1e-5 * 3.0
        flow_in_stem = sum(
            ground_area_flux.(
                z_root_depths,
                z_bottom_stem,
                p_soilf,
                Y[1],
                a_root,
                b_root,
                K_max_root,
                RAI,
            ),
        )
        flow_out_stem = ground_area_flux(
            z_bottom_stem,
            z_leaf,
            Y[1],
            Y[2],
            a_stem,
            b_stem,
            K_max_stem,
            SAI,
        )
        F[1] = flow_in_stem - flow_out_stem
        F[2] = flow_out_stem - Tf * LAI
    end


    # Check that the final state is in the new equilibrium
    soln = nlsolve(f2!, [p_leaf[end], p_stem[end]]; ftol = 1e-12)
    p_stem_f = soln.zero[1]
    p_leaf_f = soln.zero[2]
    @test abs((p_stem_f - p_stem[end]) / p_stem_f) < 1e-3
    @test abs((p_leaf_f - p_leaf[end]) / p_leaf_f) < 1e-3
end

export volumetric_liquid_fraction,
    effective_saturation,
    matric_potential,
    pressure_head,
    hydraulic_conductivity,
    impedance_factor,
    viscosity_factor

"""
    volumetric_liquid_fraction(ϑ_l::FT, ν_eff::FT) where {FT}

A pointwise function returning the volumetric liquid fraction
given the augmented liquid fraction and the effective porosity.
"""
function volumetric_liquid_fraction(ϑ_l, ν_eff::FT) where {FT}
    return min(ϑ_l, ν_eff)
end


"""
     matric_potential(α::FT, n::FT, m::FT, S::FT) where {FT}

A point-wise function returning the matric potential, using the
van Genuchten formulation.
"""
function matric_potential(α::FT, n::FT, m::FT, S) where {FT}
    ψ_m = -((S^(-FT(1) / m) - FT(1)) * α^(-n))^(FT(1) / n)
    return ψ_m
end

"""
    effective_saturation(porosity::FT, ϑ_l::FT, θr::FT) where {FT}

A point-wise function computing the effective saturation.
"""
function effective_saturation(porosity::FT, ϑ_l, θr::FT) where {FT}
    ϑ_l_safe = max(ϑ_l, θr + eps(FT))
    S_l = (ϑ_l_safe - θr) / (porosity - θr)
    return S_l
end

"""
    pressure_head(
        α::FT,
        n::FT,
        m::FT,
        θ_r::FT,
        ϑ_l::FT,
        ν_eff::FT,
        S_s::FT,
    ) where {FT}

A point-wise function returning the pressure head in
variably saturated soil, using the van Genuchten formulation
for matric potential if the soil is not saturated, and
an approximation of the positive pressure in the soil if the
soil is saturated.
"""
function pressure_head(
    α::FT,
    n::FT,
    m::FT,
    θ_r::FT,
    ϑ_l,
    ν_eff::FT,
    S_s::FT,
) where {FT}
    S_l_eff = effective_saturation(ν_eff, ϑ_l, θ_r)
    if S_l_eff <= FT(1.0)
        ψ = matric_potential(α, n, m, S_l_eff)
    else
        ψ = (ϑ_l - ν_eff) / S_s
    end
    return ψ
end

"""
     hydraulic_conductivity(K_sat::FT, m::FT, S::FT) where {FT}

A point-wise function returning the hydraulic conductivity, using the
van Genuchten formulation.
"""
function hydraulic_conductivity(K_sat, m, S::T) where {T}
    if S < 1
        K = sqrt(S) * (1 - (1 - S^(1 / m))^m)^2
    else
        K = T(1)
    end
    return K * K_sat
end


"""
    impedance_factor(
        f_i::FT,
        Ω::FT
    ) where {FT}
Returns the multiplicative factor reducing conductivity when 
a fraction of ice `f_i` is present.

Only for use with the `EnergyHydrology` model.
"""
function impedance_factor(f_i::FT, Ω::FT) where {FT}
    gamma = FT(10.0^(-Ω * f_i))
    return gamma
end


"""
    viscosity_factor(
        T::FT,
        γ::FT,
        γT_ref::FT,
    ) where {FT}

Returns the multiplicative factor which accounts for
the temperature dependence of the conductivity.

Only for use with the `EnergyHydrology` model.
"""
function viscosity_factor(T::FT, γ::FT, γT_ref::FT) where {FT}
    factor = FT(γ * (T - γT_ref))
    Theta = FT(exp(factor))
    return Theta
end

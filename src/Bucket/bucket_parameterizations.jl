"""
     heaviside(x::FT)::FT where {FT}

Computes the heaviside function.
"""
function heaviside(x::FT)::FT where {FT}
    if x > eps(FT)
        return FT(1.0)
    else
        return FT(0.0)
    end
end
"""
    melt_timescale(σS::FT, LH_f0::FT, τc::FT, F_sfc::FT, E::FT) where {FT}

Returns the timescale for melting.

This is linked to the net energy flux in the case where the surface is
warming and the melt rate is not too fast compared to the amount of snow,
and is equal to τ_c otherwise.

To be revisited.
"""
function melt_timescale(σS::FT, LH_f0::FT, τc::FT, F_sfc::FT, E::FT) where {FT}
    # if F_sfc - LH_f0*E < 0 - > net flux into land
    # so put a minus sign in τ definition
    return max(τc, -σS * LH_f0 / (F_sfc + LH_f0 * E))
end

"""
    partition_surface_fluxes(
        σS::FT,
        T_sfc::FT,
        τc::FT,
        snow_cover_fraction::FT,
        E::FT,
        F_sfc::FT,
        _LH_f0::FT,
        _T_freeze::FT,
    ) where{FT}

Partitions the surface fluxes in a flux for melting snow, a flux for sublimating snow,
and a ground heat flux. 

All fluxes are positive if they are in the direction from land towards
the atmosphere.
"""
function partition_surface_fluxes(
    σS::FT,
    T_sfc::FT,
    τc::FT,
    snow_cover_fraction::FT,
    E::FT,
    F_sfc::FT,
    _LH_f0::FT,
    _T_freeze::FT,
) where {FT}
    τ = melt_timescale(σS, _LH_f0, τc, F_sfc, E) # Eqn 23
    F_melt = -σS * _LH_f0 * heaviside(T_sfc - _T_freeze) / τ # Eqn 22. Negative by definition (into the snow/land).
    F_into_snow = -_LH_f0 * E * snow_cover_fraction + F_melt # F_melt is already multiplied by σS. Eqn 20
    G = (F_sfc - F_into_snow) # Eqn 20
    return (; F_melt = F_melt, F_into_snow = F_into_snow, G = G)
end

"""
    surface_albedo(albedo::BulkAlbedo{FT}, coords, S::FT, S_c::FT)::FT where {FT <: AbstractFloat}

Returns the bulk surface albedo, linearly interpolating between the albedo
of snow and of soil, based on the snow water equivalent S relative to
the parameter S_c.

The linear interpolation is taken from Lague et al 2019.
"""
function surface_albedo(
    albedo::BulkAlbedo{FT},
    coords,
    σS::FT,
    σS_c::FT,
)::FT where {FT <: AbstractFloat}
    (; α_snow, α_soil) = albedo
    α_soil_values = α_soil(coords)
    safe_σS::FT = max(σS, eps(FT))
    return (FT(1.0) - σS / (σS + σS_c)) * α_soil_values +
           σS / (σS + σS_c) * α_snow
end

"""
    infiltration_at_point(W::FT, M::FT, P::FT, E::FT, W_f::FT)::FT where {FT <: AbstractFloat}

Returns the soil infiltration given the current water content of the bucket W, 
the snow melt volume flux M, the precipitation volume flux P, the liquid evaporative volume
flux E, and the bucket capacity W_f. Positive values indicate increasing
soil moisture; the soil infiltration is the magnitude of the water
flux into the soil.

Extra inflow when the bucket is at capacity runs off.
Note that P and M are positive by definition, while E can be 
positive (evaporation) or negative (condensation).
"""
function infiltration_at_point(
    W::FT,
    M::FT,
    P::FT,
    E::FT,
    W_f::FT,
)::FT where {FT <: AbstractFloat}
    if (W > W_f) & (P + M - E > 0) # Equation (4b) of text
        return FT(0)
    else
        return FT(P + M - E)
    end
end

"""
    β(W::FT, W_f::FT) where {FT}

Returns the coefficient which scales the saturated
specific humidity at the surface based on the bucket water
levels, which is then used
 to obtain the
true specific humidity of the soil surface <= q_sat.
"""
function β(W::FT, W_f::FT) where {FT}
    safe_W = max(FT(0.0), W)
    if safe_W < FT(0.75) * W_f
        safe_W / (FT(0.75) * W_f)
    else
        FT(1.0)
    end
end

"""
    saturation_specific_humidity(T::FT, σS::FT, ρ_sfc::FT, parameters::PE)::FT where {FT, PE}

Computes the saturation specific humidity for the land surface, over ice
if snow is present (σS>0), and over water for a snow-free surface.
"""
function saturation_specific_humidity(
    T::FT,
    σS::FT,
    ρ_sfc::FT,
    parameters::PE,
)::FT where {FT, PE}
    thermo_params = LSMP.thermodynamic_parameters(parameters.earth_param_set)
    return (FT(1.0) - heaviside(σS)) *
           Thermodynamics.q_vap_saturation_generic(
        thermo_params,
        T,
        ρ_sfc,
        Thermodynamics.Liquid(),
    ) +
           heaviside(σS) * Thermodynamics.q_vap_saturation_generic(
        thermo_params,
        T,
        ρ_sfc,
        Thermodynamics.Ice(),
    )
end

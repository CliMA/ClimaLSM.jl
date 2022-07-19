module Parameters

abstract type AbstractLSMParameters end
const ALSMP = AbstractLSMParameters

Base.@kwdef struct LSMParameters{FT, TP, SFP} <: ALSMP
    K_therm::FT
    ρ_cloud_liq::FT
    ρ_cloud_ice::FT
    cp_l::FT
    cp_i::FT
    T_0::FT
    LH_v0::FT
    LH_s0::FT
    Stefan::FT
    T_freeze::FT
    thermo_params::TP
    surf_flux_params::SFP
end

# wrapper methods:

K_therm(ps::ALSMP) = ps.K_therm
ρ_cloud_liq(ps::ALSMP) = ps.ρ_cloud_liq
ρ_cloud_ice(ps::ALSMP) = ps.ρ_cloud_ice
cp_l(ps::ALSMP) = ps.cp_l
cp_i(ps::ALSMP) = ps.cp_i
T_0(ps::ALSMP) = ps.T_0
LH_v0(ps::ALSMP) = ps.LH_v0
LH_s0(ps::ALSMP) = ps.LH_s0
Stefan(ps::ALSMP) = ps.Stefan
T_freeze(ps::ALSMP) = ps.T_freeze
# Derived parameters
LH_f0(ps::ALSMP) = LH_s0(ps) - LH_v0(ps)

# Dependency parameter wrappers
thermodynamic_parameters(ps::ALSMP) = ps.thermo_params
surface_fluxes_parameters(ps::ALSMP) = ps.surf_flux_params

end # module

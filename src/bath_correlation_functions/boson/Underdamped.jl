export Boson_Underdamped_Matsubara

@doc raw"""
    Boson_Underdamped_Matsubara(op, λ, Γ, ω0, kT, N)
Construct an underdamped bosonic bath with Matsubara expansion

# Parameters
- `op` : The system coupling operator, must be Hermitian and, for fermionic systems, even-parity to be compatible with charge conservation.
- `λ::Real`: The coupling strength between the system and the bath.
- `Γ::Real`: The band-width of the bath spectral density.
- `ω0::Real`: The resonance frequency of the bath spectral density.
- `kT::Real`: The product of the Boltzmann constant ``k`` and the absolute temperature ``T`` of the bath.
- `N::Int`: (N+2)-terms of exponential terms are used to approximate the bath correlation function.

# Returns
- `bath::BosonBath` : a bosonic bath object with describes the interaction between system and bosonic bath
"""
function Boson_Underdamped_Matsubara(op, λ::Real, Γ::Real, ω0::Real, kT::Real, N::Int)
    Ω = sqrt(ω0^2 - (Γ/2)^2)
    ν = 2π*kT*(1:N)

    η_real = ComplexF64[(λ^2/(4*Ω))*coth((Ω + im*Γ/2)/(2*kT)), (λ^2/(4*Ω))*coth((Ω - im*Γ/2)/(2*kT))]
    γ_real = ComplexF64[Γ/2 - im*Ω, Γ/2 + im*Ω]
    η_imag = ComplexF64[(λ^2/(4*Ω))*im, -(λ^2/(4*Ω))*im]
    γ_imag = ComplexF64[Γ/2 - im*Ω, Γ/2 + im*Ω]

    if N > 0
        for l in 1:N
            append!(η_real, -2*λ^2*Γ*kT*ν[l]/(((Ω + im*Γ/2)^2 + ν[l]^2)*((Ω - im*Γ/2)^2 + ν[l]^2)))
            append!(γ_real, ν[l])
        end
    end

    return BosonBath(op, η_real, γ_real, η_imag, γ_imag)
end
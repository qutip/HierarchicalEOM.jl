"""
# `Boson_DrudeLorentz_Pade(op, λ, W, T, N)`
Constructing Drude-Lorentz bosonic bath with Padé expansion

A Padé approximant is a sum-over-poles expansion (see https://en.wikipedia.org/wiki/Pad%C3%A9_approximant).

The application of the Padé method to spectrum decompoisitions is described in Ref. [1].

[1] J. Chem. Phys. 134, 244106 (2011); https://doi.org/10.1063/1.3602466

## Parameters
- `op::AbstractMatrix` : The system operator according to the system-fermionic-bath interaction.
- `λ::Real`: The coupling strength between the system and the bath.
- `W::Real`: The reorganization energy (band-width) of the bath.
- `T::Real`: The temperature of the bath.
- `N::Int`: Number of exponential terms used to approximate the bath correlation functions.

## Returns
- `bath::BosonBath` : a bosonic bath object with describes the interaction between system and bosonic bath
"""
function Boson_DrudeLorentz_Pade(
        op::AbstractMatrix,
        λ::Real,
        W::Real,
        T::Real,
        N::Int
    )
    
    β = 1. / T
    κ, ϵ = pade_NmN(N, fermion=false)

    η = [λ * W * (cot(W * β / 2.0) - 1.0im)]
    γ = [W]
    if N > 0
        for l in 2:(N + 1)
            append!(η,
                (κ[l] / β) * 4 * λ * W * (ϵ[l] / β) / (((ϵ[l] ^ 2) / (β ^ 2)) - W ^ 2)
            )
            append!(γ, ϵ[l] / β)
        end
    end

    δ = 2 * λ / (β * W) - 1.0im * λ
    for k in 1:(N + 1)
        δ -= η[k] / γ[k]
    end

    return BosonBath(op, η, γ, δ)
end
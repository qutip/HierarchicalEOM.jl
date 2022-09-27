function f_approx(x, κ, ϵ, N)
    f = 0.5
    for l in 2:(N + 1)
        f -= 2 * κ[l] * x / (x ^ 2 + ϵ[l] ^ 2)
    end
    return f
end

function lorentzian_pade_corr(σ::Real, Γ::Real, μ::Real, W::Real, T::Real, N::Int)
    β = 1. / T
    κ, ϵ = pade_NmN(N);
    η_list = [0.5 * Γ * W * f_approx(1.0im * β * W, κ, ϵ, N)]
    γ_list = [W - σ * 1.0im * μ]

    if N > 0
        for l in 2:(N + 1)
            append!(η_list, -1.0im * (κ[l] / β) * Γ * W ^ 2 / (-(ϵ[l] / β) ^ 2 + W ^ 2))
            append!(γ_list, ϵ[l] / β - σ * 1.0im * μ)
        end
    end
    return η_list, γ_list
end

"""
# `Lorentzian_Pade_Bath(op, Γ, μ, W, T, N)`
Constructing Lorentzian fermionic bath with Padé expansion

A Padé approximant is a sum-over-poles expansion (see https://en.wikipedia.org/wiki/Pad%C3%A9_approximant).

The application of the Padé method to spectrum decompoisitions is described in Ref. [1].

[1] J. Chem. Phys. 134, 244106 (2011); https://doi.org/10.1063/1.3602466

## Parameters
- `op::AbstractMatrix` : The system operator according to the system-fermionic-bath interaction.
- `Γ::Real`: The coupling strength between the system and the bath.
- `μ::Real`: The chemical potential of the bath.
- `W::Real`: The band-width of the bath.
- `T::Real`: The temperature of the bath.
- `N::Int`: Number of exponential terms used to approximate the bath correlation functions.

## Returns
- `bath::FermionBath` : a fermionic bath object with describes the interaction between system and fermionic bath
"""
function Lorentzian_Pade_Bath(
        op::AbstractMatrix,
        Γ::Real,
        μ::Real,
        W::Real,
        T::Real,
        N::Int
    )
    η_ab, γ_ab = lorentzian_pade_corr( 1.0, Γ, μ, W, T, N)
    η_em, γ_em = lorentzian_pade_corr(-1.0, Γ, μ, W, T, N)

    return FermionBath(op, η_ab, γ_ab, η_em, γ_em)
end
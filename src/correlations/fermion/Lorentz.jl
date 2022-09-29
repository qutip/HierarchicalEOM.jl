function _fermion_lorentz_matsubara_param(σ::Real, λ::Real, μ::Real, W::Real, T::Real, N::Int)
    β = 1. / T
    ϵ = mastubara(N, fermion=true)

    η = [0.5 * λ * W * _fermi(1.0im * β * W)]
    γ = [W - σ * 1.0im * μ]

    if N > 0
        for l in 2:(N + 1)
            append!(η, 
                -1.0im * λ * T * W ^ 2 / (-(ϵ[l] * T) ^ 2 + W ^ 2)
            )
            append!(γ, ϵ[l] * T - σ * 1.0im * μ)
        end
    end
    return η, γ
end

"""
# `Fermion_Lorentz_Matsubara(op, λ, μ, W, T, N)`
Constructing Lorentzian fermionic bath with Matsubara expansion

## Parameters
- `op::AbstractMatrix` : The system \"emission\" operator according to the system-fermionic-bath interaction.
- `λ::Real`: The coupling strength between the system and the bath.
- `μ::Real`: The chemical potential of the bath.
- `W::Real`: The reorganization energy (band-width) of the bath.
- `T::Real`: The temperature of the bath.
- `N::Int`: Number of exponential terms used to approximate the bath correlation functions.

## Returns
- `bath::FermionBath` : a fermionic bath object with describes the interaction between system and fermionic bath
"""
function Fermion_Lorentz_Matsubara(
        op::AbstractMatrix,
        λ::Real,
        μ::Real,
        W::Real,
        T::Real,
        N::Int
    )
    η_ab, γ_ab = _fermion_lorentz_matsubara_param( 1.0, λ, μ, W, T, N)
    η_em, γ_em = _fermion_lorentz_matsubara_param(-1.0, λ, μ, W, T, N)

    return FermionBath(op, η_ab, γ_ab, η_em, γ_em)
end

function _fermion_lorentz_pade_param(σ::Real, λ::Real, μ::Real, W::Real, T::Real, N::Int)
    β = 1. / T
    κ, ϵ = pade_NmN(N, fermion=true)   
     
    η = [0.5 * λ * W * _fermi_pade(1.0im * β * W, κ, ϵ, N)]
    γ = [W - σ * 1.0im * μ]

    if N > 0
        for l in 2:(N + 1)
            append!(η, 
                -1.0im * κ[l] * λ * T * W ^ 2 / (-(ϵ[l] * T) ^ 2 + W ^ 2)
            )
            append!(γ, ϵ[l] * T - σ * 1.0im * μ)
        end
    end
    return η, γ
end

"""
# `Fermion_Lorentz_Pade(op, λ, μ, W, T, N)`
Constructing Lorentzian fermionic bath with Padé expansion

A Padé approximant is a sum-over-poles expansion (see https://en.wikipedia.org/wiki/Pad%C3%A9_approximant).

The application of the Padé method to spectrum decompoisitions is described in Ref. [1].

[1] J. Chem. Phys. 134, 244106 (2011); https://doi.org/10.1063/1.3602466

## Parameters
- `op::AbstractMatrix` : The system \"emission\" operator according to the system-fermionic-bath interaction.
- `λ::Real`: The coupling strength between the system and the bath.
- `μ::Real`: The chemical potential of the bath.
- `W::Real`: The reorganization energy (band-width) of the bath.
- `T::Real`: The temperature of the bath.
- `N::Int`: Number of exponential terms used to approximate the bath correlation functions.

## Returns
- `bath::FermionBath` : a fermionic bath object with describes the interaction between system and fermionic bath
"""
function Fermion_Lorentz_Pade(
        op::AbstractMatrix,
        λ::Real,
        μ::Real,
        W::Real,
        T::Real,
        N::Int
    )
    η_ab, γ_ab = _fermion_lorentz_pade_param( 1.0, λ, μ, W, T, N)
    η_em, γ_em = _fermion_lorentz_pade_param(-1.0, λ, μ, W, T, N)

    return FermionBath(op, η_ab, γ_ab, η_em, γ_em)
end
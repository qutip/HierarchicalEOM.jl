function _boson_drude_lorentz_approx_discrepancy(λ::Real, W::Real, T::Real, N::Int, η::Vector{Ti}, γ::Vector{Tj}) where {Ti, Tj <: Number}
    δ = 2 * λ * T / W - 1.0im * λ
    for k in 1:(N + 1)
        δ -= η[k] / γ[k]
    end
    return δ
end

"""
    Boson_DrudeLorentz_Matsubara(op, λ, W, T, N)
Constructing Drude-Lorentz bosonic bath with Matsubara expansion

# Parameters
- `op` : The system operator according to the system-fermionic-bath interaction.
- `λ::Real`: The coupling strength between the system and the bath.
- `W::Real`: The reorganization energy (band-width) of the bath.
- `T::Real`: The temperature of the bath.
- `N::Int`: Number of exponential terms used to approximate the bath correlation functions.

# Returns
- `bath::BosonBath` : a bosonic bath object with describes the interaction between system and bosonic bath
"""
function Boson_DrudeLorentz_Matsubara(
        op,
        λ::Real,
        W::Real,
        T::Real,
        N::Int
    )
    β = 1. / T
    ϵ = matsubara(N, fermion=false)

    η = ComplexF64[λ * W * (cot(W * β / 2.0) - 1.0im)]
    γ = ComplexF64[W]
    
    if N > 0
        for l in 2:(N + 1)
            append!(η,
                4 * λ * W * ϵ[l] * (T ^ 2) / (((ϵ[l] * T) ^ 2) - W ^ 2)
            )
            append!(γ, ϵ[l] * T)
        end
    end

    δ = _boson_drude_lorentz_approx_discrepancy(λ, W, T, N, η, γ)

    return BosonBath(op, η, γ, δ)
end

"""
    Boson_DrudeLorentz_Pade(op, λ, W, T, N)
Constructing Drude-Lorentz bosonic bath with Padé expansion

A Padé approximant is a sum-over-poles expansion (see [here](https://en.wikipedia.org/wiki/Pad%C3%A9_approximant) for more details).

The application of the Padé method to spectrum decompoisitions is described in Ref. [1].

[1] [J. Chem. Phys. 134, 244106 (2011)](https://doi.org/10.1063/1.3602466)

# Parameters
- `op` : The system operator according to the system-fermionic-bath interaction.
- `λ::Real`: The coupling strength between the system and the bath.
- `W::Real`: The reorganization energy (band-width) of the bath.
- `T::Real`: The temperature of the bath.
- `N::Int`: Number of exponential terms used to approximate the bath correlation functions.

# Returns
- `bath::BosonBath` : a bosonic bath object with describes the interaction between system and bosonic bath
"""
function Boson_DrudeLorentz_Pade(
        op,
        λ::Real,
        W::Real,
        T::Real,
        N::Int
    )
    
    β = 1. / T
    κ, ζ = pade_NmN(N, fermion=false)

    η = ComplexF64[λ * W * (cot(W * β / 2.0) - 1.0im)]
    γ = ComplexF64[W]
    if N > 0
        for l in 2:(N + 1)
            append!(η,
                κ[l] * 4 * λ * W * ζ[l] * (T ^ 2) / (((ζ[l] * T) ^ 2) - W ^ 2)
            )
            append!(γ, ζ[l] * T)
        end
    end

    δ = _boson_drude_lorentz_approx_discrepancy(λ, W, T, N, η, γ)

    return BosonBath(op, η, γ, δ)
end
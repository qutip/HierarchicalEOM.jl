abstract type AbstractEnvBath end

"""
# `BosonicBath <: AbstractEnvBath`
An object for bosonic bath correlation

## Fields
- `c_list`: the coefficient ``c_i`` in bosonic bath correlation functions (``\\sum_i c_i e^{-\\nu_i t}``).
- `ν_list`: the coefficient ``\\nu_i`` in bosonic bath correlation functions (``\\sum_i c_i e^{-\\nu_i t}``).
- `coupOP`: the operator(s) describing the coupling between system and bosonic bath.
- `N_term`: the number of exponential-expansion term of correlation function
- `N_oper`: the number of coupling operators

## Constructor
`BosonicBath(c_list, ν_list, coupOP)`
"""
mutable struct BosonicBath <: AbstractEnvBath
    c_list
    ν_list
    coupOP
    N_term::Int
    N_oper::Int

    function BosonicBath(
            c_list::Vector{Ti},
            ν_list::Vector{Tj},
            coupOP::T
        ) where {Ti, Tj, T}

        if (Ti <: Number) && (Tj <: Number) && (T <: AbstractMatrix)
            # check if the length of c_list and ν_list are valid
            N_exp_term = length(c_list)
            N_operator = 1
            if N_exp_term != length(ν_list)
                error("The length of \'c_list\' and \'ν_list\' should be the same.")
            end
        else
            error("The type of \'c_list\', \'ν_list\', and coupOP are not correct.")
        end
        return new(c_list, ν_list, coupOP, N_exp_term, N_operator)
    end
end

"""
# `FermionicBath <: AbstractEnvBath`
An object for fermionic bath correlation

## Fields
- `η_list`: the coefficient ``\\eta_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `γ_list`: the coefficient ``\\gamma_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `coupOP`: the operator(s) describing the coupling between system and fermionic bath.
- `N_term`: the number of exponential-expansion term of correlation function
- `N_oper`: the number of coupling operators

## Constructor
`FermionicBath(c_list, ν_list, coupOP)`
"""
mutable struct FermionicBath <: AbstractEnvBath
    η_list
    γ_list
    coupOP
    N_term::Int
    N_oper::Int

    function FermionicBath(
            η_list::Vector{Ti},
            γ_list::Vector{Tj},
            coupOP::T
        ) where {Ti, Tj, T}

        if (Ti <: Vector) && (Tj <: Vector) && (T <: Vector)
            if (eltype(Ti) <: Number) && (eltype(Tj) <: Number) && (T == AbstractMatrix)
                # check if the length of η_list, γ_list, and coupOP are valid
                N_operator = length(coupOP)
                N_exp_term = length(η_list[1])
                if (N_operator != length(η_list)) || (N_operator != length(γ_list))
                    error("The length of \'η_list\', \'γ_list\', and \'Coup_Ops\' should all be the same.")
                end
                for i in 1:N_operator
                    if length(η_list[i]) != N_exp_term
                        error("The length of each vector in \'η_list\' are wrong.")
                    end
                    if length(γ_list[i]) != N_exp_term
                        error("The length of each vector in \'γ_list\' are wrong.")
                    end
                end
            end
        else
            error("The type of \'η_list\', \'γ_list\', and coupOP are not correct.")
        end
        return new(η_list, γ_list, coupOP, N_exp_term, N_operator)
    end
end
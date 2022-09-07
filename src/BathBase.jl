abstract type AbstractBath end

"""
# `BosonicBath <: AbstractBath`
An object for bosonic bath correlation

## Fields
- `η_list` : the coefficient ``\\eta_i`` in bosonic bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `γ_list` : the coefficient ``\\gamma_i`` in bosonic bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `coupOP` : the operator describing the coupling between system and bosonic bath.
- `N_term` : the number of exponential-expansion term of correlation function
- `N_oper` : the number of coupling operator (should be 1 for boson)

## Constructor
`BosonicBath(η_list, γ_list, coupOP)`
"""
mutable struct BosonicBath <: AbstractBath
    η_list
    γ_list
    coupOP
    N_term::Int
    N_oper::Int

    function BosonicBath(
            η_list::Vector{Ti},
            γ_list::Vector{Tj},
            coupOP::T
        ) where {Ti, Tj, T}

        if (Ti <: Number) && (Tj <: Number) && (T <: AbstractMatrix)
            # check if the length of η_list and γ_list are valid
            N_exp_term = length(η_list)
            if N_exp_term != length(γ_list)
                error("The length of \'η_list\' and \'γ_list\' should be the same.")
            end
        else
            error("The type of \'η_list\', \'γ_list\', and coupOP are not correct.")
        end
        return new(η_list, γ_list, coupOP, N_exp_term, 1)
    end
end

"""
# `FermionicBath <: AbstractBath`
An object for fermionic bath correlation

## Fields
- `η_list` : the coefficient ``\\eta_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `γ_list` : the coefficient ``\\gamma_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `coupOP` : the operators describing the coupling between system and fermionic bath.
- `N_term` : the number of exponential-expansion term of correlation function
- `N_oper` : the number of coupling operators

## Constructor
`FermionicBath(η_list, γ_list, coupOP)`
"""
mutable struct FermionicBath <: AbstractBath
    η_list
    γ_list
    coupOP
    N_term::Int
    N_oper::Int

    function FermionicBath(
        η_list::Vector{Ti},
        γ_list::Vector{Tj},
        coupOP::Vector{Tk}
    ) where {Ti, Tj, Tk}

    if (Ti <: AbstractVector) && (eltype(Ti) <: Number) && (Tj <: AbstractVector) && (eltype(Tj) <: Number) && (Tk <: AbstractMatrix)
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
    else
        error("The type of \'η_list\', \'γ_list\', and coupOP are not correct.")
    end
    return new(η_list, γ_list, coupOP, N_exp_term, N_operator)
end
end
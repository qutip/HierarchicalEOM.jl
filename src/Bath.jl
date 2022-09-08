abstract type AbstractBath end

"""
# `CoupOp`
An object which describes the specific coupling operator between system and bath

## Fields
- `Op` : the operator acting on system which describes the coupling between system and bath.
- `dim` : the dimension of the operator (should be equal to the system dimension).
- `η_list` : the coefficient ``\\eta_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `γ_list` : the coefficient ``\\gamma_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `N_term` : the number of exponential-expansion term of correlation function

## Constructor
`CoupOp(Op, η_list, γ_list)`

- `Op::AbstractMatrix` : the operator acting on system which describes the coupling between system and bath.
- `η_list::Vector{Ti<:Number}` : the coefficient ``\\eta_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `γ_list::Vector{Tj<:Number}` : the coefficient ``\\gamma_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
"""
struct CoupOp
    Op::AbstractMatrix
    dim::Int
    η_list::AbstractVector
    γ_list::AbstractVector
    N_term::Int
    
    function CoupOp(
            Op::AbstractMatrix,
            η_list::Vector{Ti},
            γ_list::Vector{Tj}
        ) where {Ti, Tj <: Number}

        N1, N2 = size(Op)
        if N1 != N2
            error("The operator \"Op\" should be an squared matrix.")
        end

        # check if the length of η_list and γ_list are valid
        N_exp_term = length(η_list)
        if N_exp_term != length(γ_list)
            error("The length of \'η_list\' and \'γ_list\' should be the same.")
        end
        return new(Op, N1, η_list, γ_list, N_exp_term)
    end
end

"""
# `BosonBath <: AbstractBath`
An object describing the interaction between system and bosonic bath

## Fields
- `Op` : the operator acting on system which describes the coupling between system and bosonic bath.
- `dim` : the dimension of the operator (should be equal to the system dimension).
- `η_list` : the coefficient ``\\eta_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `γ_list` : the coefficient ``\\gamma_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `N_term` : the number of exponential-expansion term of correlation function

## Constructor
`BosonBath(coup::CoupOp)`
"""
struct BosonBath <: AbstractBath
    Op::AbstractMatrix
    dim::Int
    η_list::AbstractVector
    γ_list::AbstractVector
    N_term::Int
    
    function BosonBath(
            Op::AbstractMatrix,
            η_list::Vector{Ti},
            γ_list::Vector{Tj}
        ) where {Ti, Tj <: Number}

        N1, N2 = size(Op)
        if N1 != N2
            error("The operator \"Op\" should be an squared matrix.")
        end

        # check if the length of η_list and γ_list are valid
        N_exp_term = length(η_list)
        if N_exp_term != length(γ_list)
            error("The length of \'η_list\' and \'γ_list\' should be the same.")
        end
        return new(sparse(Op), N1, η_list, γ_list, N_exp_term)
    end

    function BosonBath(coup::CoupOp)
        new(sparse(coup.Op), coup.dim, coup.η_list, coup.γ_list, coup.N_term)
    end
end

"""
# `FermionBath <: AbstractBath`
An object describing the interaction between system and fermionic bath

## Fields
- `Op` : the operators acting on system which describes the coupling between system and fermionic bath.
- `dim` : the dimension of the operator (should be equal to the system dimension).
- `η_list` : the coefficient ``\\eta_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `γ_list` : the coefficient ``\\gamma_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `N_term` : the number of exponential-expansion term of correlation function

## Constructor
- `FermionBath(coup::CoupOp...)`
- `FermionBath(coup::Vector{CoupOp})`
"""
struct FermionBath <: AbstractBath
    Op::AbstractVector
    dim::Int
    η_list::AbstractVector
    γ_list::AbstractVector
    N_term::Int
    
    function FermionBath(
            Op::Vector{Ti},
            η_list::Vector{Tj},
            γ_list::Vector{Tk}
        ) where {Ti, Tj, Tk}

        if (Ti <: AbstractMatrix) && (Tj <: AbstractVector) && (eltype(Tj) <: Number) && (Tk <: AbstractVector) && (eltype(Tk) <: Number)
            # check if the length of η_list, γ_list, and Op are valid
            N, N1 = size(Op[1])
            if N != N1
                error("All the operators \"Op\" should be squared matrices with same dimension.")
            end
            N_operator = length(Op)
            N_exp_term = length(η_list[1])
            if (N_operator != length(η_list)) || (N_operator != length(γ_list))
                error("The length of \'η_list\', \'γ_list\', and \'Coup_Ops\' should all be the same.")
            end
            for i in 1:N_operator
                if (N, N) != size(Op[i])
                    error("All the operators \"Op\" should be squared matrices with same dimension.")
                end
                if length(η_list[i]) != N_exp_term
                    error("The length of each vector in \'η_list\' are wrong.")
                end
                if length(γ_list[i]) != N_exp_term
                    error("The length of each vector in \'γ_list\' are wrong.")
                end
            end
        else
            error("The type of \'η_list\', \'γ_list\', and Op are not correct.")
        end
        return new(Op, N, η_list, γ_list, N_exp_term)
    end

    function FermionBath(coup::CoupOp...) FermionBath([coup...]) end

    function FermionBath(coup::Vector{CoupOp})
        Op = [sparse(coup[1].Op)]
        dim = coup[1].dim
        η_list = [coup[1].η_list]
        γ_list = [coup[1].γ_list]
        N_term = coup[1].N_term
        
        for i in 2:length(coup)
            if coup[i].dim != dim
                error("The system dimension of the coupling operators should be equal.")
            elseif coup[i].N_term != N_term
                error("The number of exponential-expansion terms of correlation function should be equal.")
            end
            push!(Op, sparse(coup[i].Op))
            push!(η_list, coup[i].η_list)
            push!(γ_list, coup[i].γ_list)
        end
        return new(Op, dim, η_list, γ_list, N_term)
    end
end
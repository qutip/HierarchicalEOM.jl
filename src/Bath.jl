abstract type AbstractBath end

spre(q::AbstractMatrix)  = sparse(kron(Matrix(I, size(q)[1], size(q)[1]), q))
spost(q::AbstractMatrix) = sparse(kron(transpose(q), Matrix(I, size(q)[1], size(q)[1])))

"""
# `BosonBath <: AbstractBath`
An object which describes the interaction between system and bosonic bath

## Fields
- `spre`  : the super-operator (right side operator multiplication) for the coupling operator.
- `spost` : the super-operator (left side operator multiplication) for the coupling operator.
- `comm`  : the super-operator (commutator) for the coupling operator.
- `dim` : the dimension of the coupling operator (should be equal to the system dimension).
- `η` : the coefficients ``\\eta_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `γ` : the coefficients ``\\gamma_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `Nterm` : the number of exponential-expansion term of correlation function

## Constructor
`BosonBath(Op, η, γ)`

- `Op::AbstractMatrix` : The system operator according to the system-bosonic-bath interaction.
- `η::Vector{Ti<:Number}` : the coefficients ``\\eta_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `γ::Vector{Tj<:Number}` : the coefficients ``\\gamma_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
"""
struct BosonBath <: AbstractBath
    spre::AbstractMatrix
    spost::AbstractMatrix
    comm::AbstractMatrix
    dim::Int
    η::AbstractVector
    γ::AbstractVector
    Nterm::Int
    
    function BosonBath(
            Op::AbstractMatrix,
            η::Vector{Ti},
            γ::Vector{Tj}
        ) where {Ti, Tj <: Number}

        N1, N2 = size(Op)
        if N1 != N2
            error("The operator \"Op\" should be an squared matrix.")
        end

        # check if the length of coefficients are valid
        N_exp_term = length(η)
        if N_exp_term != length(γ)
            error("The length of \'η\' and \'γ\' should be the same.")
        end
        spreQ  = spre(Op)
        spostQ = spost(Op)
        return new(spreQ, spostQ, spreQ - spostQ, N1, η, γ, N_exp_term)
    end
end

"""
# `FermionBath <: AbstractBath`
An object which describes the interaction between system and fermionic bath

## Fields
- `spre`   : the super-operator (right side operator multiplication) for the coupling operator.
- `spost`  : the super-operator (left side operator multiplication) for the coupling operator.
- `spreD`  : the super-operator (right side operator multiplication) for the adjoint of the coupling operator.
- `spostD` : the super-operator (left side operator multiplication) for the adjoint of the coupling operator.
- `dim` : the dimension of the coupling operator (should be equal to the system dimension).
- `η_absorb` : the coefficients ``\\eta_i`` for absorption in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `γ_absorb` : the coefficients ``\\gamma_i`` for absorption in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `η_emit` : the coefficients ``\\eta_i`` for emission in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `γ_emit` : the coefficients ``\\gamma_i`` for emission in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `Nterm` : the number of exponential-expansion term of correlation function

## Constructor
`FermionBath(Op, η_list, γ_list)`

- `Op::AbstractMatrix` : The system operator according to the system-fermionic-bath interaction.
- `η_absorb` : the coefficients ``\\eta_i`` for absorption in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `γ_absorb` : the coefficients ``\\gamma_i`` for absorption in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `η_emit` : the coefficients ``\\eta_i`` for emission in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `γ_emit` : the coefficients ``\\gamma_i`` for emission in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
"""
struct FermionBath <: AbstractBath
    spre::AbstractMatrix
    spost::AbstractMatrix
    spreD::AbstractMatrix
    spostD::AbstractMatrix
    dim::Int
    η_absorb::AbstractVector
    γ_absorb::AbstractVector
    η_emit::AbstractVector
    γ_emit::AbstractVector
    Nterm::Int
    
    function FermionBath(
            Op::AbstractMatrix,
            η_absorb::Vector{Tk},
            γ_absorb::Vector{Tl},
            η_emit::Vector{Ti},
            γ_emit::Vector{Tj}
        ) where {Ti, Tj, Tk, Tl <: Number}

        N1, N2 = size(Op)
        if N1 != N2
            error("The operator \"Op\" should be an squared matrix.")
        end

        # check if the length of coefficients are valid
        N_exp_term = length(η_emit)
        if (N_exp_term != length(γ_emit)) || (N_exp_term != length(η_absorb)) || (N_exp_term != length(γ_absorb))
            error("The length of \'η_emit\', \'γ_emit\', \'η_absorb\' and \'γ_absorb\' should all be the same.")
        end
        return new(spre(Op), spost(Op), spre(adjoint(Op)), spost(adjoint(Op)), N1, η_absorb, γ_absorb, η_emit, γ_emit, N_exp_term)
    end
end
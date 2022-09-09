abstract type AbstractBath end

spre(q::AbstractMatrix)  = sparse(kron(Matrix(I, size(q)[1], size(q)[1]), q))
spost(q::AbstractMatrix) = sparse(kron(transpose(q), Matrix(I, size(q)[1], size(q)[1])))

"""
# `bosonOP`
An object which describes the coupling operator between system and bosonic bath

## Fields
- `spre`  : the super-operator (right side operator multiplication) for the operator.
- `spost` : the super-operator (left side operator multiplication) for the operator.
- `comm`  : the super-operator (commutator) for the operator.
- `dim` : the dimension of the operator (should be equal to the system dimension).
- `η_list` : the coefficient ``\\eta_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `γ_list` : the coefficient ``\\gamma_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `N_term` : the number of exponential-expansion term of correlation function

## Constructor
`bosonOp(Op, η_list, γ_list)`

- `Op::AbstractMatrix` : The system operator according to the system-bosonic-bath interaction.
- `η_list::Vector{Ti<:Number}` : the coefficient ``\\eta_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `γ_list::Vector{Tj<:Number}` : the coefficient ``\\gamma_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
"""
struct bosonOP
    spre::AbstractMatrix
    spost::AbstractMatrix
    comm::AbstractMatrix
    dim::Int
    η_list::AbstractVector
    γ_list::AbstractVector
    N_term::Int
    
    function bosonOP(
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
        spreQ  = spre(Op)
        spostQ = spost(Op)
        return new(spreQ, spostQ, spreQ - spostQ, N1, η_list, γ_list, N_exp_term)
    end
end

"""
# `BosonBath <: AbstractBath`
An object describing the interaction between system and bosonic bath

## Fields
- `coup::bosonOP` : the object which describes the coupling operator between system and bosonic bath
- `dim` : the dimension of the system dimension.
- `N_term` : the total number of exponential-expansion term of correlation function

## Constructor
`BosonBath(coup::bosonOP)`
"""
struct BosonBath <: AbstractBath
    coup::bosonOP
    dim::Int
    N_term::Int
    
    function BosonBath(
            Op::AbstractMatrix,
            η_list::Vector{Ti},
            γ_list::Vector{Tj}
        ) where {Ti, Tj <: Number}

        bop = bosonOP(Op, η_list, γ_list)
        return new(bop, bop.dim, bop.N_term)
    end

    function BosonBath(coup::bosonOP)
        new(coup, coup.dim, coup.N_term)
    end
end

"""
# `fermionOP`
An object which describes the coupling operator between system and fermionic bath

## Fields
- `spre`   : the super-operator (right side operator multiplication) for the operator.
- `spost`  : the super-operator (left side operator multiplication) for the operator.
- `spreD`  : the super-operator (right side operator multiplication) for the adjoint of the operator.
- `spostD` : the super-operator (left side operator multiplication) for the adjoint of the operator.
- `dim` : the dimension of the operator (should be equal to the system dimension).
- `η_list` : the coefficient ``\\eta_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `γ_list` : the coefficient ``\\gamma_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `N_term` : the number of exponential-expansion term of correlation function

## Constructor
`fermionOP(Op, η_list, γ_list)`

- `Op::AbstractMatrix` : The system operator according to the system-fermionic-bath interaction.
- `η_list::Vector{Ti<:Number}` : the coefficient ``\\eta_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `γ_list::Vector{Tj<:Number}` : the coefficient ``\\gamma_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
"""
struct fermionOP
    spre::AbstractMatrix
    spost::AbstractMatrix
    spreD::AbstractMatrix
    sposD::AbstractMatrix
    dim::Int
    η_list::AbstractVector
    γ_list::AbstractVector
    N_term::Int
    
    function fermionOP(
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
        return new(spre(Op), spost(Op), spre(adjoint(Op)), spost(adjoint(Op)), N1, η_list, γ_list, N_exp_term)
    end
end

"""
# `FermionBath <: AbstractBath`
An object describing the interaction between system and fermionic bath

## Fields
- `coup::Vector{fermionOP}` : the operators acting on system which describes the coupling between system and fermionic bath.
- `dim` : the dimension of system.
- `N_term` : the total number of exponential-expansion term of correlation function

## Constructor
- `FermionBath(coup::fermionOP...)`
- `FermionBath(coup::Vector{fermionOP})`
"""
struct FermionBath <: AbstractBath
    coup::Vector{fermionOP}
    dim::Int
    N_term::Int
    
    function FermionBath(
            Op::Vector{Ti},
            η_list::Vector{Tj},
            γ_list::Vector{Tk}
        ) where {Ti, Tj, Tk}

        if (Ti <: AbstractMatrix) && (Tj <: AbstractVector) && (eltype(Tj) <: Number) && (Tk <: AbstractVector) && (eltype(Tk) <: Number)
            dim, = size(Op[1])
            N_operator = length(Op)
            if (N_operator != length(η_list)) || (N_operator != length(γ_list))
                error("The length of \'η_list\', \'γ_list\', and \'Coup_Ops\' should all be the same.")
            end

            Coup = []
            N_exp_term = 0
            for i in 1:N_operator
                coup = fermionOP(Op[i], η_list[i], γ_list[i])
                if coup.dim != dim
                    error("All the dimensions of operators should be the same.")
                end
                push!(Coup, coup)
                N_exp_term += coup.N_term
            end
        else
            error("The type of \'η_list\', \'γ_list\', and \'Op\' are not correct.")
        end
        return new(Coup, dim, N_exp_term)
    end

    function FermionBath(coup::fermionOP...) FermionBath([coup...]) end

    function FermionBath(coup::Vector{fermionOP})
        dim = coup[1].dim
        N_exp_term = 0
        N_operator = length(coup)
        for i in 1:N_operator
            if coup[i].dim != dim
                error("All the dimensions of operators should be the same.")
            end
            N_exp_term += coup.N_term
        end

        return new(coup, dim, N_exp_term)
    end
end
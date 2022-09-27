abstract type AbstractBath end
abstract type AbstractBosonBath end
abstract type AbstractFermionBath end

spre(q::AbstractMatrix)  = sparse(kron(Matrix(I, size(q)[1], size(q)[1]), q))
spost(q::AbstractMatrix) = sparse(kron(transpose(q), Matrix(I, size(q)[1], size(q)[1])))

show(io::IO, B::AbstractBath) = print(io, "$(typeof(B)) object with (system) dim = $(B.dim) and $(B.Nterm) exponential-expansion terms\n")
show(io::IO, m::MIME"text/plain", B::AbstractBath) =  show(io, B)

show(io::IO, B::AbstractBosonBath) = print(io, "$(typeof(B))-type bath with (system) dim = $(B.dim) and $(B.Nterm) exponential-expansion terms\n")
show(io::IO, m::MIME"text/plain", B::AbstractBosonBath) =  show(io, B)

show(io::IO, B::AbstractFermionBath) = print(io, "$(typeof(B))-type bath with (system) dim = $(B.dim) and $(B.Nterm) exponential-expansion terms\n")
show(io::IO, m::MIME"text/plain", B::AbstractFermionBath) =  show(io, B)

isclose(a::Number, b::Number, rtol=1e-05, atol=1e-08) = abs(a - b) <= (atol + rtol * abs(b))

function combine_same_gamma(η::Vector{Ti}, γ::Vector{Tj}) where {Ti, Tj <: Number}
    if length(η) != length(γ)
        error("The length of \'η\' and \'γ\' should be the same.")
    end

    ηnew  = ComplexF64[η[1]]
    γnew  = ComplexF64[γ[1]]
    for j in 2:length(γ)
        valueDontExist = true
        for k in 1:length(γnew)
            if isclose(γnew[k], γ[j])
                ηnew[k] += η[j]
                valueDontExist = false
            end
        end
        if valueDontExist
            push!(ηnew, η[j])
            push!(γnew, γ[j])
        end
    end
    return ηnew, γnew
end

"""
# `BosonBath <: AbstractBath`
An object which describes the interaction between system and bosonic bath

## Fields
- `bath::Vector{T<:AbstractBosonBath}` : the vector for different operators which describes the interaction between system and bosonic bath
- `dim` : the dimension of the coupling operator (should be equal to the system dimension).
- `Nterm` : the number of exponential-expansion term of correlation functions

## Constructor
### `BosonBath(op, η, γ, combine)`
For the case where real part and imaginary part of the correlation function are combined.

- `op::AbstractMatrix` : The system operator according to the system-bosonic-bath interaction.
- `η::Vector{Ti<:Number}` : the coefficients ``\\eta_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `γ::Vector{Tj<:Number}` : the coefficients ``\\gamma_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `combine::Bool` : Whether to combine the exponential-expansion terms with the same frequency. Defaults to `true`.

### `BosonBath(op, η_real, γ_real, η_imag, γ_imag, combine)`
For the case where the correlation function splits into real part and imaginary part.

- `op::AbstractMatrix` : The system operator according to the system-bosonic-bath interaction.
- `η_real::Vector{Ti<:Number}` : the real part of coefficients ``\\eta_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `γ_real::Vector{Tj<:Number}` : the real part of coefficients ``\\gamma_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `η_imag::Vector{Tk<:Number}` : the imaginary part of coefficients ``\\eta_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `γ_imag::Vector{Tl<:Number}` : the imaginary part of coefficients ``\\gamma_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `combine::Bool` : Whether to combine the exponential-expansion terms with the same frequency. Defaults to `true`.
"""
struct BosonBath <: AbstractBath
    bath::AbstractVector
    dim::Int
    Nterm::Int

    function BosonBath(baths::Vector{T}) where T <: AbstractBosonBath
        dim   = baths[1].dim
        Nterm = baths[1].Nterm
        if length(baths) > 1
            for b in baths[2:end]
                if b.dim != dim 
                    error("The system dimension of the bosonic bath coupling operators are not consistent.")
                end
                Nterm += b.Nterm
            end
        else
            return new(baths, dim, Nterm)
        end
        return new(baths, dim, Nterm)
    end

    function BosonBath(
            op::AbstractMatrix,
            η::Vector{Ti},
            γ::Vector{Tj},
            combine=true
        ) where {Ti, Tj <: Number}
        if combine
            ηnew, γnew = combine_same_gamma(η, γ)
            bRI = bosonRealImag(op, real.(ηnew), imag.(ηnew), γnew)
        else
            bRI = bosonRealImag(op, real.(η), imag.(η), γ)
        end
        return new([bRI], bRI.dim, bRI.Nterm)
    end

    function BosonBath(
            op::AbstractMatrix,
            η_real::Vector{Ti},
            γ_real::Vector{Tj},
            η_imag::Vector{Tk},
            γ_imag::Vector{Tl},
            combine::Bool=true
        ) where {Ti, Tj, Tk, Tl <: Number}

        if combine
            ηR, γR = combine_same_gamma(η_real, γ_real)
            ηI, γI = combine_same_gamma(η_imag, γ_imag)

            NR = length(γR)
            NI = length(γI)
            Nterm = NR + NI

            ηRI_real = ComplexF64[]
            ηRI_imag = ComplexF64[]
            γRI      = ComplexF64[]
            real_idx = Int64[]
            imag_idx = Int64[]
            for j in 1:NR
                for k in 1:NI
                    if isclose(γR[j],  γI[k])
                        # record index
                        push!(real_idx, j)
                        push!(imag_idx, k)

                        # append combined coefficient vectors
                        push!(ηRI_real, ηR[j])
                        push!(ηRI_imag, ηI[k])
                        push!(γRI, γR[j])
                    end
                end
            end
            sort!(real_idx)
            sort!(imag_idx)
            deleteat!(ηR, real_idx)
            deleteat!(γR, real_idx)
            deleteat!(ηI, imag_idx)
            deleteat!(γI, imag_idx)
            
            bR  = bosonReal(op, ηR, γR)
            bI  = bosonImag(op, ηI, γI)
            bRI = bosonRealImag(op, ηRI_real, ηRI_imag, γRI)
            Nterm_new = bR.Nterm + bI.Nterm + bRI.Nterm
            if Nterm != (Nterm_new + bRI.Nterm)
                error("Conflicts occur in combining real and imaginary parts of bath correlation function.")
            end
            return new([bR, bI, bRI], bR.dim, Nterm_new)

        else
            bR = bosonReal(op, η_real, γ_real)
            bI = bosonImag(op, η_imag, γ_imag)
            return new([bR, bI], bR.dim, bR.Nterm + bI.Nterm)
        end
    end
end

combineBath(bath::BosonBath...) = combineBath([bath...])

function combineBath(bath::Vector{BosonBath})
    baths = AbstractBosonBath[]
    for b in bath
        push!(baths, b.bath...)
    end
    return BosonBath(baths)
end

"""
# `bosonReal <: AbstractBosonBath`
A bosonic bath for the real part of bath correlation function

## Fields
- `Comm`  : the super-operator (commutator) for the coupling operator.
- `dim` : the dimension of the coupling operator (should be equal to the system dimension).
- `η` : the real part of coefficients ``\\eta_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `γ` : the real part of coefficients ``\\gamma_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `Nterm` : the number of exponential-expansion term of correlation function

## Constructor
`bosonReal(op, η_real, γ_real)`

- `op::AbstractMatrix` : The system operator according to the system-bosonic-bath interaction.
- `η_real::Vector{Ti<:Number}` : the real part of coefficients ``\\eta_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `γ_real::Vector{Tj<:Number}` : the real part of coefficients ``\\gamma_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
"""
struct bosonReal <: AbstractBosonBath
    Comm::AbstractMatrix
    dim::Int
    η::AbstractVector
    γ::AbstractVector
    Nterm::Int
    
    function bosonReal(
            op::AbstractMatrix,
            η_real::Vector{Ti},
            γ_real::Vector{Tj}
        ) where {Ti, Tj <: Number}

        N1, N2 = size(op)
        if N1 != N2
            error("The operator \"op\" should be an squared matrix.")
        end

        # check if the length of coefficients are valid
        N_exp_term = length(η_real)
        if N_exp_term != length(γ_real)
            error("The length of \'η_real\' and \'γ_real\' should be the same.")
        end

        return new(spre(op) - spost(op), N1, η_real, γ_real, N_exp_term)
    end
end

"""
# `bosonImag <: AbstractBosonBath`
A bosonic bath for the imaginary part of bath correlation function

## Fields
- `Comm`  : the super-operator (commutator) for the coupling operator.
- `anComm`  : the super-operator (anti-commutator) for the coupling operator.
- `dim` : the dimension of the coupling operator (should be equal to the system dimension).
- `η` : the imaginary part of coefficients ``\\eta_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `γ` : the imaginary part of coefficients ``\\gamma_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `Nterm` : the number of exponential-expansion term of correlation function

## Constructor
`bosonImag(op, η_imag, γ_imag)`

- `op::AbstractMatrix` : The system operator according to the system-bosonic-bath interaction.
- `η_imag::Vector{Ti<:Number}` : the imaginary part of coefficients ``\\eta_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `γ_imag::Vector{Tj<:Number}` : the imaginary part of coefficients ``\\gamma_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
"""
struct bosonImag <: AbstractBosonBath
    Comm::AbstractMatrix
    anComm::AbstractMatrix
    dim::Int
    η::AbstractVector
    γ::AbstractVector
    Nterm::Int
    
    function bosonImag(
            op::AbstractMatrix,
            η_imag::Vector{Ti},
            γ_imag::Vector{Tj}
        ) where {Ti, Tj <: Number}

        N1, N2 = size(op)
        if N1 != N2
            error("The operator \"op\" should be an squared matrix.")
        end

        # check if the length of coefficients are valid
        N_exp_term = length(η_imag)
        if N_exp_term != length(γ_imag)
            error("The length of \'η_imag\' and \'γ_imag\' should be the same.")
        end
        spreQ  = spre(op)
        spostQ = spost(op)
        return new(spreQ - spostQ, spreQ + spostQ, N1, η_imag, γ_imag, N_exp_term)
    end
end

"""
# `bosonRealImag <: AbstractBosonBath`
A bosonic bath which the real part and imaginary part of the bath correlation function are combined

## Fields
- `Comm`  : the super-operator (commutator) for the coupling operator.
- `anComm`  : the super-operator (anti-commutator) for the coupling operator.
- `dim` : the dimension of the coupling operator (should be equal to the system dimension).
- `η_real` : the real part of coefficients ``\\eta_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `η_imag` : the imaginary part of coefficients ``\\eta_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `γ` : the coefficients ``\\gamma_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `Nterm` : the number of exponential-expansion term of correlation function

## Constructor
`bosonRealImag(op, η_real, η_imag, γ)`

- `op::AbstractMatrix` : The system operator according to the system-bosonic-bath interaction.
- `η_real` : the real part of coefficients ``\\eta_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `η_imag` : the imaginary part of coefficients ``\\eta_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `γ::Vector{Tj<:Number}` : the coefficients ``\\gamma_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
"""
struct bosonRealImag <: AbstractBosonBath
    Comm::AbstractMatrix
    anComm::AbstractMatrix
    dim::Int
    η_real::AbstractVector
    η_imag::AbstractVector
    γ::AbstractVector
    Nterm::Int
    
    function bosonRealImag(
            op::AbstractMatrix,
            η_real::Vector{Ti},
            η_imag::Vector{Tj},
            γ::Vector{Tk}
        ) where {Ti, Tj, Tk <: Number}

        N1, N2 = size(op)
        if N1 != N2
            error("The operator \"op\" should be an squared matrix.")
        end

        # check if the length of coefficients are valid
        N_exp_term = length(η_real)
        if (N_exp_term != length(η_imag)) || (N_exp_term != length(γ))
            error("The length of \'η_real\', \'η_imag\' and \'γ\' should be the same.")
        end
        spreQ  = spre(op)
        spostQ = spost(op)
        return new(spreQ - spostQ, spreQ + spostQ, N1, η_real, η_imag, γ, N_exp_term)
    end
end

"""
# `FermionBath <: AbstractBath`
An object which describes the interaction between system and fermionic bath

## Fields
- `bath::Vector{T<:AbstractFermionBath}` : the vector for different operators which describes the interaction
- `dim` : the dimension of the coupling operator (should be equal to the system dimension).
- `Nterm` : the number of exponential-expansion term of correlation functions

## Constructor
### `FermionBath(op, η_absorb, γ_absorb, η_emit, γ_emit)`

- `op::AbstractMatrix` : The system \"emission\" operator according to the system-fermionic-bath interaction.
- `η_absorb` : the coefficients ``\\eta_i`` for absorption in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `γ_absorb` : the coefficients ``\\gamma_i`` for absorption in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `η_emit` : the coefficients ``\\eta_i`` for emission in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `γ_emit` : the coefficients ``\\gamma_i`` for emission in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
"""
struct FermionBath <: AbstractBath
    bath::AbstractVector
    dim::Int
    Nterm::Int

    function FermionBath(baths::Vector{T}) where T <: AbstractFermionBath
        dim   = baths[1].dim
        Nterm = baths[1].Nterm
        if length(baths) > 1
            for b in baths[2:end]
                if b.dim != dim 
                    error("The system dimension of the fermionic bath coupling operators are not consistent.")
                end
                Nterm += b.Nterm
            end
        else
            return new(baths, dim, Nterm)
        end
        return new(baths, dim, Nterm)
    end 

    function FermionBath(
            op::AbstractMatrix,
            η_absorb::Vector{Ti},
            γ_absorb::Vector{Tj},
            η_emit::Vector{Tk},
            γ_emit::Vector{Tl}
        ) where {Ti, Tj, Tk, Tl <: Number}

        fA = fermionAbsorb(adjoint(op), η_absorb, γ_absorb, η_emit)
        fE = fermionEmit(op, η_emit, γ_emit, η_absorb)
        return new([fA, fE], fA.dim, fA.Nterm + fE.Nterm)
    end
end

combineBath(bath::FermionBath...) = combineBath([bath...])

function combineBath(bath::Vector{FermionBath})
    baths = AbstractFermionBath[]
    for b in bath
        push!(baths, b.bath...)
    end
    return FermionBath(baths)
end

"""
# `fermionAbsorb <: AbstractFermionBath`
An object which describes the absorption of fermion in the interaction

## Fields
- `spre`   : the super-operator (right side operator multiplication) for the coupling operator.
- `spost`  : the super-operator (left side operator multiplication) for the coupling operator.
- `spreD`  : the super-operator (right side operator multiplication) for the adjoint of the coupling operator.
- `spostD` : the super-operator (left side operator multiplication) for the adjoint of the coupling operator.
- `dim` : the dimension of the coupling operator (should be equal to the system dimension).
- `η` : the coefficients ``\\eta_i`` for absorption in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `γ` : the coefficients ``\\gamma_i`` for absorption in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `η_emit` : the coefficients ``\\eta_i`` for emission in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `Nterm` : the number of exponential-expansion term of correlation function

## Constructor
`fermionAbsorb(op, η_absorb, γ_absorb, η_emit)`

- `op::AbstractMatrix` : The system absorption operator according to the system-fermionic-bath interaction.
- `η_absorb` : the coefficients ``\\eta_i`` for absorption in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `γ_absorb` : the coefficients ``\\gamma_i`` for absorption in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `η_emit` : the coefficients ``\\eta_i`` for emission in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
"""
struct fermionAbsorb <: AbstractFermionBath
    spre::AbstractMatrix
    spost::AbstractMatrix
    spreD::AbstractMatrix
    spostD::AbstractMatrix
    dim::Int
    η::AbstractVector
    γ::AbstractVector
    η_emit::AbstractVector
    Nterm::Int
    
    function fermionAbsorb(
            op::AbstractMatrix,
            η_absorb::Vector{Ti},
            γ_absorb::Vector{Tj},
            η_emit::Vector{Tk}
        ) where {Ti, Tj, Tk <: Number}

        N1, N2 = size(op)
        if N1 != N2
            error("The operator \"op\" should be an squared matrix.")
        end

        # check if the length of coefficients are valid
        N_exp_term = length(η_absorb)
        if (N_exp_term != length(γ_absorb)) || (N_exp_term != length(η_emit))
            error("The length of \'η_absorb\', \'γ_absorb\' and \'η_emit\' should all be the same.")
        end
        return new(spre(op), spost(op), spre(adjoint(op)), spost(adjoint(op)), N1, η_absorb, γ_absorb, η_emit, N_exp_term)
    end
end

"""
# `fermionEmit <: AbstractFermionBath`
An object which describes the emission of fermion in the interaction

## Fields
- `spre`   : the super-operator (right side operator multiplication) for the coupling operator.
- `spost`  : the super-operator (left side operator multiplication) for the coupling operator.
- `spreD`  : the super-operator (right side operator multiplication) for the adjoint of the coupling operator.
- `spostD` : the super-operator (left side operator multiplication) for the adjoint of the coupling operator.
- `dim` : the dimension of the coupling operator (should be equal to the system dimension).
- `η` : the coefficients ``\\eta_i`` for emission in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `γ` : the coefficients ``\\gamma_i`` for emission in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `η_absorb` : the coefficients ``\\eta_i`` for absorption in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `Nterm` : the number of exponential-expansion term of correlation function

## Constructor
`fermionEmit(op, η_emit, γ_emit, η_absorb)`

- `op::AbstractMatrix` : The system emission operator according to the system-fermionic-bath interaction.
- `η_emit` : the coefficients ``\\eta_i`` for emission in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `γ_emit` : the coefficients ``\\gamma_i`` for emission in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `η_absorb` : the coefficients ``\\eta_i`` for absorption in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
"""
struct fermionEmit <: AbstractFermionBath
    spre::AbstractMatrix
    spost::AbstractMatrix
    spreD::AbstractMatrix
    spostD::AbstractMatrix
    dim::Int
    η::AbstractVector
    γ::AbstractVector
    η_absorb::AbstractVector
    Nterm::Int
    
    function fermionEmit(
            op::AbstractMatrix,
            η_emit::Vector{Ti},
            γ_emit::Vector{Tj},
            η_absorb::Vector{Tk}
        ) where {Ti, Tj, Tk <: Number}

        N1, N2 = size(op)
        if N1 != N2
            error("The operator \"op\" should be an squared matrix.")
        end

        # check if the length of coefficients are valid
        N_exp_term = length(η_emit)
        if (N_exp_term != length(γ_emit)) || (N_exp_term != length(η_absorb))
            error("The length of \'η_emit\', \'γ_emit\' and \'η_absorb\' should all be the same.")
        end
        return new(spre(op), spost(op), spre(adjoint(op)), spost(adjoint(op)), N1, η_emit, γ_emit, η_absorb, N_exp_term)
    end
end
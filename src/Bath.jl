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

function _check_coupling_operator(op)
    if isValidMatrixType(op)
        N1, N2 = size(op)
        if N1 != N2
            error("The operator \"op\" should be an squared matrix.")
        end
        return N1
    else
        error("Invalid matrix type of op.")
    end
end

function _combine_same_gamma(η::Vector{Ti}, γ::Vector{Tj}) where {Ti, Tj <: Number}
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
    struct BosonBath <: AbstractBath
An object which describes the interaction between system and bosonic bath

# Fields
- `bath` : the different boson-bath-type objects which describes the interaction between system and bosonic bath
- `op` : The system operator according to the system-bosonic-bath interaction.
- `dim` : the dimension of the coupling operator (should be equal to the system dimension).
- `Nterm` : the number of exponential-expansion term of correlation functions
- `δ` : The approximation discrepancy which is used for adding the terminator to HEOM matrix (see function: addTerminator!)
"""
struct BosonBath <: AbstractBath
    bath::Vector{AbstractBosonBath}
    op::AbstractMatrix
    dim::Int
    Nterm::Int
    δ::Number
end

"""
    BosonBath(op, η, γ, δ=0.0; combine=true)
Generate BosonBath object for the case where real part and imaginary part of the correlation function are combined.

# Parameters
- `op` : The system operator according to the system-bosonic-bath interaction.
- `η::Vector{Ti<:Number}` : the coefficients ``\\eta_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `γ::Vector{Tj<:Number}` : the coefficients ``\\gamma_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `δ::Number` : The approximation discrepancy (Default to `0.0`) which is used for adding the terminator to HEOM matrix (see function: addTerminator!)
- `combine::Bool` : Whether to combine the exponential-expansion terms with the same frequency. Defaults to `true`.
"""
function BosonBath(
        op,
        η::Vector{Ti},
        γ::Vector{Tj},
        δ::Number=0.0;
        combine::Bool=true
    ) where {Ti, Tj <: Number}
    if combine
        ηnew, γnew = _combine_same_gamma(η, γ)
        bRI = bosonRealImag(op, real.(ηnew), imag.(ηnew), γnew)
    else
        bRI = bosonRealImag(op, real.(η), imag.(η), γ)
    end
    return BosonBath(AbstractBosonBath[bRI], copy(op), bRI.dim, bRI.Nterm, δ)
end

"""
    BosonBath(op, η_real, γ_real, η_imag, γ_imag, δ=0.0; combine=true)
Generate BosonBath object for the case where the correlation function splits into real part and imaginary part.

# Parameters
- `op` : The system operator according to the system-bosonic-bath interaction.
- `η_real::Vector{Ti<:Number}` : the real part of coefficients ``\\eta_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `γ_real::Vector{Tj<:Number}` : the real part of coefficients ``\\gamma_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `η_imag::Vector{Tk<:Number}` : the imaginary part of coefficients ``\\eta_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `γ_imag::Vector{Tl<:Number}` : the imaginary part of coefficients ``\\gamma_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `δ::Number` : The approximation discrepancy (Default to `0.0`) which is used for adding the terminator to HEOM matrix (see function: addTerminator!)
- `combine::Bool` : Whether to combine the exponential-expansion terms with the same frequency. Defaults to `true`.
"""
function BosonBath(
        op,
        η_real::Vector{Ti},
        γ_real::Vector{Tj},
        η_imag::Vector{Tk},
        γ_imag::Vector{Tl},
        δ::Tm=0.0;
        combine::Bool=true
    ) where {Ti, Tj, Tk, Tl, Tm <: Number}

    if combine
        ηR, γR = _combine_same_gamma(η_real, γ_real)
        ηI, γI = _combine_same_gamma(η_imag, γ_imag)

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
        return BosonBath(AbstractBosonBath[bR, bI, bRI], copy(op), bR.dim, Nterm_new, δ)

    else
        bR = bosonReal(op, η_real, γ_real)
        bI = bosonImag(op, η_imag, γ_imag)
        return BosonBath(AbstractBosonBath[bR, bI], copy(op), bR.dim, bR.Nterm + bI.Nterm, δ)
    end
end

"""
    struct bosonReal <: AbstractBosonBath
A bosonic bath for the real part of bath correlation function

# Fields
- `Comm`  : the super-operator (commutator) for the coupling operator.
- `dim` : the dimension of the coupling operator (should be equal to the system dimension).
- `η` : the real part of coefficients ``\\eta_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `γ` : the real part of coefficients ``\\gamma_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `Nterm` : the number of exponential-expansion term of correlation function
"""
struct bosonReal <: AbstractBosonBath
    Comm::SparseMatrixCSC{ComplexF64, Int64}
    dim::Int
    η::AbstractVector
    γ::AbstractVector
    Nterm::Int
end

"""
    bosonReal(op, η_real, γ_real)
Generate bosonic bath for the real part of bath correlation function

# Parameters
- `op` : The system operator according to the system-bosonic-bath interaction.
- `η_real::Vector{Ti<:Number}` : the real part of coefficients ``\\eta_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `γ_real::Vector{Tj<:Number}` : the real part of coefficients ``\\gamma_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
"""
function bosonReal(
        op,
        η_real::Vector{Ti},
        γ_real::Vector{Tj}
    ) where {Ti, Tj <: Number}

    dim = _check_coupling_operator(op)

    # check if the length of coefficients are valid
    N_exp_term = length(η_real)
    if N_exp_term != length(γ_real)
        error("The length of \'η_real\' and \'γ_real\' should be the same.")
    end

    return bosonReal(spre(op) - spost(op), dim, η_real, γ_real, N_exp_term)
end

"""
    struct bosonImag <: AbstractBosonBath
A bosonic bath for the imaginary part of bath correlation function

# Fields
- `Comm`  : the super-operator (commutator) for the coupling operator.
- `anComm`  : the super-operator (anti-commutator) for the coupling operator.
- `dim` : the dimension of the coupling operator (should be equal to the system dimension).
- `η` : the imaginary part of coefficients ``\\eta_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `γ` : the imaginary part of coefficients ``\\gamma_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `Nterm` : the number of exponential-expansion term of correlation function
"""
struct bosonImag <: AbstractBosonBath
    Comm::SparseMatrixCSC{ComplexF64, Int64}
    anComm::SparseMatrixCSC{ComplexF64, Int64}
    dim::Int
    η::AbstractVector
    γ::AbstractVector
    Nterm::Int
end

"""
    bosonImag(op, η_imag, γ_imag)
Generate bosonic bath for the imaginary part of correlation function

# Parameters
- `op` : The system operator according to the system-bosonic-bath interaction.
- `η_imag::Vector{Ti<:Number}` : the imaginary part of coefficients ``\\eta_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `γ_imag::Vector{Tj<:Number}` : the imaginary part of coefficients ``\\gamma_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
"""
    function bosonImag(
        op,
        η_imag::Vector{Ti},
        γ_imag::Vector{Tj}
    ) where {Ti, Tj <: Number}

    dim = _check_coupling_operator(op)

    # check if the length of coefficients are valid
    N_exp_term = length(η_imag)
    if N_exp_term != length(γ_imag)
        error("The length of \'η_imag\' and \'γ_imag\' should be the same.")
    end
    spreQ  = spre(op)
    spostQ = spost(op)
    return bosonImag(spreQ - spostQ, spreQ + spostQ, dim, η_imag, γ_imag, N_exp_term)
end

"""
    sturct bosonRealImag <: AbstractBosonBath
A bosonic bath which the real part and imaginary part of the bath correlation function are combined

# Fields
- `Comm`  : the super-operator (commutator) for the coupling operator.
- `anComm`  : the super-operator (anti-commutator) for the coupling operator.
- `dim` : the dimension of the coupling operator (should be equal to the system dimension).
- `η_real` : the real part of coefficients ``\\eta_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `η_imag` : the imaginary part of coefficients ``\\eta_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `γ` : the coefficients ``\\gamma_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `Nterm` : the number of exponential-expansion term of correlation function
"""
struct bosonRealImag <: AbstractBosonBath
    Comm::SparseMatrixCSC{ComplexF64, Int64}
    anComm::SparseMatrixCSC{ComplexF64, Int64}
    dim::Int
    η_real::AbstractVector
    η_imag::AbstractVector
    γ::AbstractVector
    Nterm::Int
end

"""
    bosonRealImag(op, η_real, η_imag, γ)
Generate bosonic bath which the real part and imaginary part of the bath correlation function are combined

# Parameters
- `op` : The system operator according to the system-bosonic-bath interaction.
- `η_real::Vector{Ti<:Number}` : the real part of coefficients ``\\eta_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `η_imag::Vector{Tj<:Number}` : the imaginary part of coefficients ``\\eta_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `γ::Vector{Tk<:Number}` : the coefficients ``\\gamma_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
"""
function bosonRealImag(
        op,
        η_real::Vector{Ti},
        η_imag::Vector{Tj},
        γ::Vector{Tk}
    ) where {Ti, Tj, Tk <: Number}

    dim = _check_coupling_operator(op)

    # check if the length of coefficients are valid
    N_exp_term = length(η_real)
    if (N_exp_term != length(η_imag)) || (N_exp_term != length(γ))
        error("The length of \'η_real\', \'η_imag\' and \'γ\' should be the same.")
    end
    spreQ  = spre(op)
    spostQ = spost(op)
    return bosonRealImag(spreQ - spostQ, spreQ + spostQ, dim, η_real, η_imag, γ, N_exp_term)
end

"""
    struct FermionBath <: AbstractBath
An object which describes the interaction between system and fermionic bath

# Fields
- `bath` : the different fermion-bath-type objects which describes the interaction
- `op` : The system \"emission\" operator according to the system-fermionic-bath interaction.
- `dim` : the dimension of the coupling operator (should be equal to the system dimension).
- `Nterm` : the number of exponential-expansion term of correlation functions
- `δ` : The approximation discrepancy which is used for adding the terminator to HEOM matrix (see function: addTerminator!)
"""
struct FermionBath <: AbstractBath
    bath::Vector{AbstractFermionBath}
    op::AbstractMatrix
    dim::Int
    Nterm::Int
    δ::Number
end

"""
    FermionBath(op, η_absorb, γ_absorb, η_emit, γ_emit, δ=0.0)
Generate FermionBath object

# Parameters
- `op` : The system \"emission\" operator according to the system-fermionic-bath interaction.
- `η_absorb::Vector{Ti<:Number}` : the coefficients ``\\eta_i`` for absorption in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `γ_absorb::Vector{Tj<:Number}` : the coefficients ``\\gamma_i`` for absorption in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `η_emit::Vector{Tk<:Number}` : the coefficients ``\\eta_i`` for emission in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `γ_emit::Vector{Tl<:Number}` : the coefficients ``\\gamma_i`` for emission in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `δ::Number` : The approximation discrepancy (Defaults to `0.0`) which is used for adding the terminator to HEOM matrix (see function: addTerminator!)
"""
function FermionBath(
        op,
        η_absorb::Vector{Ti},
        γ_absorb::Vector{Tj},
        η_emit::Vector{Tk},
        γ_emit::Vector{Tl},
        δ::Tm=0.0
    ) where {Ti, Tj, Tk, Tl, Tm <: Number}

    fA = fermionAbsorb(adjoint(op), η_absorb, γ_absorb, η_emit)
    fE = fermionEmit(op, η_emit, γ_emit, η_absorb)
    return FermionBath(AbstractFermionBath[fA, fE], copy(op), fA.dim, fA.Nterm + fE.Nterm, δ)
end

"""
    struct fermionAbsorb <: AbstractFermionBath
An object which describes the absorption of the system in the interaction

# Fields
- `spre`   : the super-operator (right side operator multiplication) for the coupling operator.
- `spost`  : the super-operator (left side operator multiplication) for the coupling operator.
- `spreD`  : the super-operator (right side operator multiplication) for the adjoint of the coupling operator.
- `spostD` : the super-operator (left side operator multiplication) for the adjoint of the coupling operator.
- `dim` : the dimension of the coupling operator (should be equal to the system dimension).
- `η` : the coefficients ``\\eta_i`` for absorption in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `γ` : the coefficients ``\\gamma_i`` for absorption in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `η_emit` : the coefficients ``\\eta_i`` for emission in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `Nterm` : the number of exponential-expansion term of correlation function
"""
struct fermionAbsorb <: AbstractFermionBath
    spre::SparseMatrixCSC{ComplexF64, Int64}
    spost::SparseMatrixCSC{ComplexF64, Int64}
    spreD::SparseMatrixCSC{ComplexF64, Int64}
    spostD::SparseMatrixCSC{ComplexF64, Int64}
    dim::Int
    η::AbstractVector
    γ::AbstractVector
    η_emit::AbstractVector
    Nterm::Int
end

"""
    fermionAbsorb(op, η_absorb, γ_absorb, η_emit)
Generate fermionic bath which describes the absorption of the system in the interaction

# Parameters
- `op` : The system absorption operator according to the system-fermionic-bath interaction.
- `η_absorb::Vector{Ti<:Number}` : the coefficients ``\\eta_i`` for absorption in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `γ_absorb::Vector{Tj<:Number}` : the coefficients ``\\gamma_i`` for absorption in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `η_emit::Vector{Tk<:Number}` : the coefficients ``\\eta_i`` for emission in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
"""
function fermionAbsorb(
        op,
        η_absorb::Vector{Ti},
        γ_absorb::Vector{Tj},
        η_emit::Vector{Tk}
    ) where {Ti, Tj, Tk <: Number}

    dim = _check_coupling_operator(op)

    # check if the length of coefficients are valid
    N_exp_term = length(η_absorb)
    if (N_exp_term != length(γ_absorb)) || (N_exp_term != length(η_emit))
        error("The length of \'η_absorb\', \'γ_absorb\' and \'η_emit\' should all be the same.")
    end
    return fermionAbsorb(spre(op), spost(op), spre(adjoint(op)), spost(adjoint(op)), dim, η_absorb, γ_absorb, η_emit, N_exp_term)
end

"""
    struct fermionEmit <: AbstractFermionBath
An object which describes the emission of the system in the interaction

# Fields
- `spre`   : the super-operator (right side operator multiplication) for the coupling operator.
- `spost`  : the super-operator (left side operator multiplication) for the coupling operator.
- `spreD`  : the super-operator (right side operator multiplication) for the adjoint of the coupling operator.
- `spostD` : the super-operator (left side operator multiplication) for the adjoint of the coupling operator.
- `dim` : the dimension of the coupling operator (should be equal to the system dimension).
- `η` : the coefficients ``\\eta_i`` for emission in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `γ` : the coefficients ``\\gamma_i`` for emission in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `η_absorb` : the coefficients ``\\eta_i`` for absorption in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `Nterm` : the number of exponential-expansion term of correlation function
"""
struct fermionEmit <: AbstractFermionBath
    spre::SparseMatrixCSC{ComplexF64, Int64}
    spost::SparseMatrixCSC{ComplexF64, Int64}
    spreD::SparseMatrixCSC{ComplexF64, Int64}
    spostD::SparseMatrixCSC{ComplexF64, Int64}
    dim::Int
    η::AbstractVector
    γ::AbstractVector
    η_absorb::AbstractVector
    Nterm::Int
end

"""
    fermionEmit(op, η_emit, γ_emit, η_absorb)
Generate fermionic bath which describes the absorption of the system in the interaction

# Parameters
- `op` : The system emission operator according to the system-fermionic-bath interaction.
- `η_emit::Vector{Ti<:Number}` : the coefficients ``\\eta_i`` for emission in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `γ_emit::Vector{Ti<:Number}` : the coefficients ``\\gamma_i`` for emission in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `η_absorb::Vector{Ti<:Number}` : the coefficients ``\\eta_i`` for absorption in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
"""
function fermionEmit(
        op,
        η_emit::Vector{Ti},
        γ_emit::Vector{Tj},
        η_absorb::Vector{Tk}
    ) where {Ti, Tj, Tk <: Number}

    dim = _check_coupling_operator(op)

    # check if the length of coefficients are valid
    N_exp_term = length(η_emit)
    if (N_exp_term != length(γ_emit)) || (N_exp_term != length(η_absorb))
        error("The length of \'η_emit\', \'γ_emit\' and \'η_absorb\' should all be the same.")
    end
    return fermionEmit(spre(op), spost(op), spre(adjoint(op)), spost(adjoint(op)), dim, η_emit, γ_emit, η_absorb, N_exp_term)
end

struct CombinedBath <: AbstractBath
    bath::AbstractVector
    dim::Int
    Nterm::Int
end

function CombinedBath(dim::Int, bath::BosonBath...)   CombinedBath(dim, [bath...]) end
function CombinedBath(dim::Int, bath::FermionBath...) CombinedBath(dim, [bath...]) end

function CombinedBath(dim::Int, B::Vector{BosonBath})
    Nterm = 0
    baths = AbstractBosonBath[]
    for b in B
        if b.dim != dim 
            error("The system dimension of the bosonic bath coupling operators are not consistent.")
        end
        push!(baths, b.bath...)
        Nterm += b.Nterm
    end

    return CombinedBath(baths, dim, Nterm)
end

function CombinedBath(dim::Int, B::Vector{FermionBath})
    Nterm = 0
    baths = AbstractFermionBath[]
    for b in B
        if b.dim != dim 
            error("The system dimension of the fermionic bath coupling operators are not consistent.")
        end
        push!(baths, b.bath...)
        Nterm += b.Nterm
    end

    return CombinedBath(baths, dim, Nterm)
end
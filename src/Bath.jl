abstract type AbstractBath end
abstract type AbstractBosonBath end
abstract type AbstractFermionBath end

@doc raw"""
    struct Exponent
An object which describes a single exponential-expansion term (naively, an excitation mode) within the decomposition of the bath correlation functions.

The expansion of a bath correlation function can be expressed as : ``C(t) = \sum_i \eta_i \exp(-\gamma_i t)``.

# Fields
- `op` : The system coupling operator according to system-bath interaction.
- `η::Number` : the coefficient ``\eta_i`` in bath correlation function.
- `γ::Number` : the coefficient ``\gamma_i`` in bath correlation function.
- `types::String` : The type-tag of the exponent.

The different types of the Exponent:
- `\"bR\"` : from real part of bosonic correlation function ``C^{u=\textrm{R}}(t)``
- `\"bI\"` : from imaginary part of bosonic correlation function ``C^{u=\textrm{I}}(t)``
- `\"bRI\"` : from combined (real and imaginary part) bosonic bath correlation function ``C(t)``
- `\"fA\"` : from absorption fermionic correlation function ``C^{\nu=+}(t)``
- `\"fE\"` : from emission fermionic correlation function ``C^{\nu=-}(t)``
"""
struct Exponent 
    op
    η::Number
    γ::Number
    types::String
end

spre(q::AbstractMatrix)  = sparse(kron(Matrix(I, size(q)[1], size(q)[1]), q))
spost(q::AbstractMatrix) = sparse(kron(transpose(q), Matrix(I, size(q)[1], size(q)[1])))

function show(io::IO, E::Exponent)
    print(io, 
        "Bath Exponent with types = \"$(E.types)\", operator size = $(size(E.op)), η = $(E.η), γ = $(E.γ).\n"
    )
end
show(io::IO, m::MIME"text/plain", E::Exponent) =  show(io, E)

show(io::IO, B::AbstractBath) = print(io, "$(typeof(B)) object with (system) dim = $(B.dim) and $(B.Nterm) exponential-expansion terms\n")
show(io::IO, m::MIME"text/plain", B::AbstractBath) =  show(io, B)

show(io::IO, B::AbstractBosonBath) = print(io, "$(typeof(B))-type bath with (system) dim = $(B.dim) and $(B.Nterm) exponential-expansion terms\n")
show(io::IO, m::MIME"text/plain", B::AbstractBosonBath) =  show(io, B)

show(io::IO, B::AbstractFermionBath) = print(io, "$(typeof(B))-type bath with (system) dim = $(B.dim) and $(B.Nterm) exponential-expansion terms\n")
show(io::IO, m::MIME"text/plain", B::AbstractFermionBath) =  show(io, B)

function checkbounds(B::AbstractBath, i::Int)
    if (i < 1) || (i > B.Nterm)
        error("Attempt to access $(B.Nterm)-exponent term Bath at index [$(i)]")
    end
end

length(B::AbstractBath) = B.Nterm
lastindex(B::AbstractBath) = B.Nterm

function getindex(B::AbstractBath, i::Int)
    checkbounds(B, i)
    
    count = 0
    for b in B.bath
        if i <= (count + b.Nterm)
            k = i - count
            b_type = typeof(b)
            op = copy(B.op)
            if b_type == bosonRealImag 
                η = b.η_real[k] + 1.0im * b.η_imag[k]
                types = "bRI"
            else
                η = b.η[k]
                if b_type == bosonReal
                    types = "bR"
                elseif b_type == bosonImag
                    types = "bI"
                elseif b_type == fermionAbsorb
                    types = "fA"
                    op = op'
                elseif b_type == fermionEmit
                    types = "fE"
                end
            end
            return Exponent(op, η, b.γ[k], types)
        else 
            count += b.Nterm
        end
    end
end

function getindex(B::AbstractBath, r::UnitRange{Int})
    checkbounds(B, r[1])
    checkbounds(B, r[end])

    count = 0
    exp_list = Exponent[]
    for b in B.bath
        for k in 1:b.Nterm
            count += 1
            if (r[1] <= count) && (count <= r[end])
                b_type = typeof(b)
                op = copy(B.op)
                if b_type == bosonRealImag 
                    η = b.η_real[k] + 1.0im * b.η_imag[k]
                    types = "bRI"
                else
                    η = b.η[k]
                    if b_type == bosonReal
                        types = "bR"
                    elseif b_type == bosonImag
                        types = "bI"
                    elseif b_type == fermionAbsorb
                        types = "fA"
                        op = op'
                    elseif b_type == fermionEmit
                        types = "fE"
                    end
                end
                push!(exp_list, Exponent(op, η, b.γ[k], types))
            end
            if count == r[end] 
                return exp_list
            end
        end
    end
end

function getindex(B::AbstractBath, ::Colon)
    return getindex(B, 1:B.Nterm)
end

iterate(B::AbstractBath) = B[1], 2
iterate(B::AbstractBath, ::Nothing) = nothing
function iterate(B::AbstractBath, state) 
    if state < length(B)
        return B[state], state + 1
    else
        return B[state], nothing
    end
end

isclose(a::Number, b::Number, rtol=1e-05, atol=1e-08) = abs(a - b) <= (atol + rtol * abs(b))

function _check_bosonic_coupling_operator(op)
    if isValidMatrixType(op)
        N,  = size(op)
        if !ishermitian(op)
            @warn "The system-bosonic-bath coupling operator \"op\" should be Hermitian operator."
        end
        return N
    else
        error("Invalid matrix \"op\".")
    end
end

function _check_fermionic_coupling_operator(op)
    if isValidMatrixType(op)
        N,  = size(op)
        return N
    else
        error("Invalid matrix \"op\".")
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

@doc raw"""
    struct BosonBath <: AbstractBath
An object which describes the interaction between system and bosonic bath

# Fields
- `bath` : the different boson-bath-type objects which describes the interaction between system and bosonic bath
- `op` : The system coupling operator, must be Hermitian and, for fermionic systems, even-parity to be compatible with charge conservation.
- `dim` : the dimension of the coupling operator (should be equal to the system dimension).
- `Nterm` : the number of exponential-expansion term of correlation functions
- `δ` : The approximation discrepancy which is used for adding the terminator to HEOM matrix (see function: addTerminator)

# Methods
One can obtain the ``k``-th exponent (exponential-expansion term) from `bath::BosonBath` by calling : `bath[k]`.
`HEOM.jl` also supports the following calls (methods) :
```julia
bath[1:k];   # returns a vector which contains the exponents from the `1`-st to the `k`-th term.
bath[1:end]; # returns a vector which contains all the exponential-expansion terms
bath[:];     # returns a vector which contains all the exponential-expansion terms
from b in bath
    # do something
end
```
"""
struct BosonBath <: AbstractBath
    bath::Vector{AbstractBosonBath}
    op::AbstractMatrix
    dim::Int
    Nterm::Int
    δ::Number
end

@doc raw"""
    BosonBath(op, η, γ, δ=0.0; combine=true)
Generate BosonBath object for the case where real part and imaginary part of the correlation function are combined.

```math
\begin{aligned}
C(\tau)
&=\frac{1}{2\pi}\int_{0}^{\infty} d\omega J(\omega)\left[n(\omega)e^{i\omega \tau}+(n(\omega)+1)e^{-i\omega \tau}\right]\\
&=\sum_i \eta_i \exp(-\gamma_i \tau),
\end{aligned}
```
where ``J(\omega)`` is the spectral density of the bath and ``n(\omega)`` represents the Bose-Einstein distribution.

# Parameters
- `op` : The system coupling operator, must be Hermitian and, for fermionic systems, even-parity to be compatible with charge conservation.
- `η::Vector{Ti<:Number}` : the coefficients ``\eta_i`` in bath correlation function ``C(\tau)``.
- `γ::Vector{Tj<:Number}` : the coefficients ``\gamma_i`` in bath correlation function ``C(\tau)``.
- `δ::Number` : The approximation discrepancy (Default to `0.0`) which is used for adding the terminator to HEOM matrix (see function: addTerminator)
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

@doc raw"""
    BosonBath(op, η_real, γ_real, η_imag, γ_imag, δ=0.0; combine=true)
Generate BosonBath object for the case where the correlation function splits into real part and imaginary part.

```math
\begin{aligned}
C(\tau)
&=\frac{1}{2\pi}\int_{0}^{\infty} d\omega J(\omega)\left[n(\omega)e^{i\omega \tau}+(n(\omega)+1)e^{-i\omega \tau}\right]\\
&=\sum_i \eta_i \exp(-\gamma_i \tau),
\end{aligned}
```
where ``J(\omega)`` is the spectral density of the bath and ``n(\omega)`` represents the Bose-Einstein distribution.

When ``\gamma_i \neq \gamma_i^*``, a closed form for the HEOM can be obtained by further decomposing ``C(\tau)`` into its real (R) and
imaginary (I) parts as
```math
C(\tau)=\sum_{u=\textrm{R},\textrm{I}}(\delta_{u, \textrm{R}} + i\delta_{u, \textrm{I}})C^{u}(\tau)
```
where ``\delta`` is the Kronecker delta function and ``C^{u}(\tau)=\sum_i \eta_i^u \exp(-\gamma_i^u \tau)``

# Parameters
- `op` : The system coupling operator, must be Hermitian and, for fermionic systems, even-parity to be compatible with charge conservation.
- `η_real::Vector{Ti<:Number}` : the coefficients ``\eta_i`` in real part of bath correlation function ``C^{u=\textrm{R}}``.
- `γ_real::Vector{Tj<:Number}` : the coefficients ``\gamma_i`` in real part of bath correlation function ``C^{u=\textrm{R}}``.
- `η_imag::Vector{Tk<:Number}` : the coefficients ``\eta_i`` in imaginary part of bath correlation function ``C^{u=\textrm{I}}``.
- `γ_imag::Vector{Tl<:Number}` : the coefficients ``\gamma_i`` in imaginary part of bath correlation function ``C^{u=\textrm{I}}``.
- `δ::Number` : The approximation discrepancy (Default to `0.0`) which is used for adding the terminator to HEOM matrix (see function: addTerminator)
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

@doc raw"""
    struct bosonReal <: AbstractBosonBath
A bosonic bath for the real part of bath correlation function ``C^{u=\textrm{R}}``

# Fields
- `Comm`  : the super-operator (commutator) for the coupling operator.
- `dim` : the dimension of the coupling operator (should be equal to the system dimension).
- `η` : the coefficients ``\eta_i`` in real part of bath correlation function ``C^{u=\textrm{R}}``.
- `γ` : the coefficients ``\gamma_i`` in real part of bath correlation function ``C^{u=\textrm{R}}``.
- `Nterm` : the number of exponential-expansion term of correlation function
"""
struct bosonReal <: AbstractBosonBath
    Comm::SparseMatrixCSC{ComplexF64, Int64}
    dim::Int
    η::AbstractVector
    γ::AbstractVector
    Nterm::Int
end

@doc raw"""
    bosonReal(op, η_real, γ_real)
Generate bosonic bath for the real part of bath correlation function ``C^{u=\textrm{R}}``

# Parameters
- `op` : The system coupling operator, must be Hermitian and, for fermionic systems, even-parity to be compatible with charge conservation.
- `η_real::Vector{Ti<:Number}` : the coefficients ``\eta_i`` in real part of bath correlation function ``C^{u=\textrm{R}}``.
- `γ_real::Vector{Tj<:Number}` : the coefficients ``\gamma_i`` in real part of bath correlation function ``C^{u=\textrm{R}}``.
"""
function bosonReal(
        op,
        η_real::Vector{Ti},
        γ_real::Vector{Tj}
    ) where {Ti, Tj <: Number}

    dim = _check_bosonic_coupling_operator(op)

    # check if the length of coefficients are valid
    N_exp_term = length(η_real)
    if N_exp_term != length(γ_real)
        error("The length of \'η_real\' and \'γ_real\' should be the same.")
    end

    return bosonReal(spre(op) - spost(op), dim, η_real, γ_real, N_exp_term)
end

@doc raw"""
    struct bosonImag <: AbstractBosonBath
A bosonic bath for the imaginary part of bath correlation function ``C^{u=\textrm{I}}``

# Fields
- `Comm`  : the super-operator (commutator) for the coupling operator.
- `anComm`  : the super-operator (anti-commutator) for the coupling operator.
- `dim` : the dimension of the coupling operator (should be equal to the system dimension).
- `η` : the coefficients ``\eta_i`` in imaginary part of bath correlation function ``C^{u=\textrm{I}}``.
- `γ` : the coefficients ``\gamma_i`` in imaginary part of bath correlation function ``C^{u=\textrm{I}}``.
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

@doc raw"""
    bosonImag(op, η_imag, γ_imag)
Generate bosonic bath for the imaginary part of correlation function ``C^{u=\textrm{I}}``

# Parameters
- `op` : The system coupling operator, must be Hermitian and, for fermionic systems, even-parity to be compatible with charge conservation.
- `η_imag::Vector{Ti<:Number}` : the coefficients ``\eta_i`` in imaginary part of bath correlation functions ``C^{u=\textrm{I}}``.
- `γ_imag::Vector{Tj<:Number}` : the coefficients ``\gamma_i`` in imaginary part of bath correlation functions ``C^{u=\textrm{I}}``.
"""
    function bosonImag(
        op,
        η_imag::Vector{Ti},
        γ_imag::Vector{Tj}
    ) where {Ti, Tj <: Number}

    dim = _check_bosonic_coupling_operator(op)

    # check if the length of coefficients are valid
    N_exp_term = length(η_imag)
    if N_exp_term != length(γ_imag)
        error("The length of \'η_imag\' and \'γ_imag\' should be the same.")
    end
    spreQ  = spre(op)
    spostQ = spost(op)
    return bosonImag(spreQ - spostQ, spreQ + spostQ, dim, η_imag, γ_imag, N_exp_term)
end

@doc raw"""
    sturct bosonRealImag <: AbstractBosonBath
A bosonic bath which the real part and imaginary part of the bath correlation function are combined 

# Fields
- `Comm`  : the super-operator (commutator) for the coupling operator.
- `anComm`  : the super-operator (anti-commutator) for the coupling operator.
- `dim` : the dimension of the coupling operator (should be equal to the system dimension).
- `η_real` : the real part of coefficients ``\eta_i`` in bath correlation function ``\sum_i \eta_i \exp(-\gamma_i t)``.
- `η_imag` : the imaginary part of coefficients ``\eta_i`` in bath correlation function ``\sum_i \eta_i \exp(-\gamma_i t)``.
- `γ` : the coefficients ``\gamma_i`` in bath correlation function ``\sum_i \eta_i \exp(-\gamma_i t)``.
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

@doc raw"""
    bosonRealImag(op, η_real, η_imag, γ)
Generate bosonic bath which the real part and imaginary part of the bath correlation function are combined

# Parameters
- `op` : The system coupling operator, must be Hermitian and, for fermionic systems, even-parity to be compatible with charge conservation.
- `η_real::Vector{Ti<:Number}` : the real part of coefficients ``\eta_i`` in bath correlation function ``\sum_i \eta_i \exp(-\gamma_i t)``.
- `η_imag::Vector{Tj<:Number}` : the imaginary part of coefficients ``\eta_i`` in bath correlation function ``\sum_i \eta_i \exp(-\gamma_i t)``.
- `γ::Vector{Tk<:Number}` : the coefficients ``\gamma_i`` in bath correlation function ``\sum_i \eta_i \exp(-\gamma_i t)``.
"""
function bosonRealImag(
        op,
        η_real::Vector{Ti},
        η_imag::Vector{Tj},
        γ::Vector{Tk}
    ) where {Ti, Tj, Tk <: Number}

    dim = _check_bosonic_coupling_operator(op)

    # check if the length of coefficients are valid
    N_exp_term = length(η_real)
    if (N_exp_term != length(η_imag)) || (N_exp_term != length(γ))
        error("The length of \'η_real\', \'η_imag\' and \'γ\' should be the same.")
    end
    spreQ  = spre(op)
    spostQ = spost(op)
    return bosonRealImag(spreQ - spostQ, spreQ + spostQ, dim, η_real, η_imag, γ, N_exp_term)
end

@doc raw"""
    struct FermionBath <: AbstractBath
An object which describes the interaction between system and fermionic bath

# Fields
- `bath` : the different fermion-bath-type objects which describes the interaction
- `op` : The system \"emission\" operator according to the system-fermionic-bath interaction.
- `dim` : the dimension of the coupling operator (should be equal to the system dimension).
- `Nterm` : the number of exponential-expansion term of correlation functions
- `δ` : The approximation discrepancy which is used for adding the terminator to HEOM matrix (see function: addTerminator)

# Methods
One can obtain the ``k``-th exponent (exponential-expansion term) from `bath::FermionBath` by calling : `bath[k]`.
`HEOM.jl` also supports the following calls (methods) :
```julia
bath[1:k];   # returns a vector which contains the exponents from the `1`-st to the `k`-th term.
bath[1:end]; # returns a vector which contains all the exponential-expansion terms
bath[:];     # returns a vector which contains all the exponential-expansion terms
from b in bath
    # do something
end
```
"""
struct FermionBath <: AbstractBath
    bath::Vector{AbstractFermionBath}
    op::AbstractMatrix
    dim::Int
    Nterm::Int
    δ::Number
end

@doc raw"""
    FermionBath(op, η_absorb, γ_absorb, η_emit, γ_emit, δ=0.0)
Generate FermionBath object

```math
\begin{aligned}
C^{\nu=+}(\tau)
&=\frac{1}{2\pi}\int_{-\infty}^{\infty} d\omega J(\omega) n(\omega) e^{i\omega \tau}\\
&=\sum_i \eta_i^{\nu=+} \exp(-\gamma_i^{\nu=+} \tau),\\
C^{\nu=-}(\tau)
&=\frac{1}{2\pi}\int_{-\infty}^{\infty} d\omega J(\omega) (1-n(\omega)) e^{-i\omega \tau}\\
&=\sum_i \eta_i^{\nu=-} \exp(-\gamma_i^{\nu=-} \tau),
\end{aligned}
```
where ``\nu=+`` (``\nu=-``) represents absorption (emission) process, ``J(\omega)`` is the spectral density of the bath and ``n(\omega)`` is the Fermi-Dirac distribution.

# Parameters
- `op` : The system annihilation operator according to the system-fermionic-bath interaction.
- `η_absorb::Vector{Ti<:Number}` : the coefficients ``\eta_i`` of absorption bath correlation function ``C^{\nu=+}(\tau)``.
- `γ_absorb::Vector{Tj<:Number}` : the coefficients ``\gamma_i`` of absorption bath correlation function ``C^{\nu=+}(\tau)``.
- `η_emit::Vector{Tk<:Number}` : the coefficients ``\eta_i`` of emission bath correlation function ``C^{\nu=-}(\tau)``.
- `γ_emit::Vector{Tl<:Number}` : the coefficients ``\gamma_i`` of emission bath correlation function ``C^{\nu=-}(\tau)``.
- `δ::Number` : The approximation discrepancy (Defaults to `0.0`) which is used for adding the terminator to HEOM matrix (see function: addTerminator)
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

@doc raw"""
    struct fermionAbsorb <: AbstractFermionBath
An bath object which describes the absorption process of the fermionic system by a correlation function ``C^{\nu=+}``

# Fields
- `spre`   : the super-operator (right side operator multiplication) for the coupling operator.
- `spost`  : the super-operator (left side operator multiplication) for the coupling operator.
- `spreD`  : the super-operator (right side operator multiplication) for the adjoint of the coupling operator.
- `spostD` : the super-operator (left side operator multiplication) for the adjoint of the coupling operator.
- `dim` : the dimension of the coupling operator (should be equal to the system dimension).
- `η` : the coefficients ``\eta_i`` of absorption bath correlation function ``C^{\nu=+}``.
- `γ` : the coefficients ``\gamma_i`` of absorption bath correlation function ``C^{\nu=+}``.
- `η_emit` : the coefficients ``\eta_i`` of emission bath correlation function ``C^{\nu=-}``.
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

@doc raw"""
    fermionAbsorb(op, η_absorb, γ_absorb, η_emit)
Generate fermionic bath which describes the absorption process of the fermionic system by a correlation function ``C^{\nu=+}``

    # Parameters
- `op` : The system creation operator according to the system-fermionic-bath interaction.
- `η_absorb::Vector{Ti<:Number}` : the coefficients ``\eta_i`` of absorption bath correlation function ``C^{\nu=+}``.
- `γ_absorb::Vector{Tj<:Number}` : the coefficients ``\gamma_i`` of absorption bath correlation function ``C^{\nu=+}``.
- `η_emit::Vector{Tk<:Number}` : the coefficients ``\eta_i`` of emission bath correlation function ``C^{\nu=-}``.
"""
function fermionAbsorb(
        op,
        η_absorb::Vector{Ti},
        γ_absorb::Vector{Tj},
        η_emit::Vector{Tk}
    ) where {Ti, Tj, Tk <: Number}

    dim = _check_fermionic_coupling_operator(op)

    # check if the length of coefficients are valid
    N_exp_term = length(η_absorb)
    if (N_exp_term != length(γ_absorb)) || (N_exp_term != length(η_emit))
        error("The length of \'η_absorb\', \'γ_absorb\' and \'η_emit\' should all be the same.")
    end
    return fermionAbsorb(spre(op), spost(op), spre(adjoint(op)), spost(adjoint(op)), dim, η_absorb, γ_absorb, η_emit, N_exp_term)
end

@doc raw"""
    struct fermionEmit <: AbstractFermionBath
An bath object which describes the emission process of the fermionic system by a correlation function ``C^{\nu=-}``

# Fields
- `spre`   : the super-operator (right side operator multiplication) for the coupling operator.
- `spost`  : the super-operator (left side operator multiplication) for the coupling operator.
- `spreD`  : the super-operator (right side operator multiplication) for the adjoint of the coupling operator.
- `spostD` : the super-operator (left side operator multiplication) for the adjoint of the coupling operator.
- `dim` : the dimension of the coupling operator (should be equal to the system dimension).
- `η` : the coefficients ``\eta_i`` of emission bath correlation function ``C^{\nu=-}``.
- `γ` : the coefficients ``\gamma_i`` of emission bath correlation function ``C^{\nu=-}``.
- `η_absorb` : the coefficients ``\eta_i`` of absorption bath correlation function ``C^{\nu=+}``.
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

@doc raw"""
    fermionEmit(op, η_emit, γ_emit, η_absorb)
Generate fermionic bath which describes the emission process of the fermionic system by a correlation function ``C^{\nu=-}``

# Parameters
- `op` : The system annihilation operator according to the system-fermionic-bath interaction.
- `η_emit::Vector{Ti<:Number}` : the coefficients ``\eta_i`` of emission bath correlation function ``C^{\nu=-}``.
- `γ_emit::Vector{Ti<:Number}` : the coefficients ``\gamma_i`` of emission bath correlation function ``C^{\nu=-}``.
- `η_absorb::Vector{Ti<:Number}` : the coefficients ``\eta_i`` of absorption bath correlation function ``C^{\nu=+}``.
"""
function fermionEmit(
        op,
        η_emit::Vector{Ti},
        γ_emit::Vector{Tj},
        η_absorb::Vector{Tk}
    ) where {Ti, Tj, Tk <: Number}

    dim = _check_fermionic_coupling_operator(op)

    # check if the length of coefficients are valid
    N_exp_term = length(η_emit)
    if (N_exp_term != length(γ_emit)) || (N_exp_term != length(η_absorb))
        error("The length of \'η_emit\', \'γ_emit\' and \'η_absorb\' should all be the same.")
    end
    return fermionEmit(spre(op), spost(op), spre(adjoint(op)), spost(adjoint(op)), dim, η_emit, γ_emit, η_absorb, N_exp_term)
end

@doc raw"""
    C(bath, tlist)
Calculate the correlation function ``C(t)`` for a given bosonic bath and time list.

Note that
```math
C(t)=\sum_{u=\textrm{R},\textrm{I}}(\delta_{u, \textrm{R}} + i\delta_{u, \textrm{I}})C^{u}(t)
```

# Parameters
- `bath::BosonBath` : The bath object which describes a certain bosonic bath.
- `tlist::AbstractVector`: The specific time.

# Returns
- `clist::Vector{ComplexF64}` : a list of the value of correlation function according to the given time list.
"""
function C(bath::BosonBath, tlist::AbstractVector)
    clist = zeros(ComplexF64, length(tlist))
    for (i, t) in enumerate(tlist)
        for e in bath
            if e.types == "bI"
                clist[i] += 1.0im * e.η * exp(- e.γ * t)
            else
                clist[i] += e.η * exp(- e.γ * t)
            end
        end
    end
    return clist
end

@doc raw"""
    C(bath, tlist)
Calculate the correlation function ``C^{\nu=+}(t)`` and ``C^{\nu=-}(t)`` for a given fermionic bath and time list.
Here, ``\nu=+`` represents the absorption process and ``\nu=-`` represents the emmision process.

# Parameters
- `bath::FermionBath` : The bath object which describes a certain fermionic bath.
- `tlist::AbstractVector`: The specific time.

# Returns
- `cplist::Vector{ComplexF64}` : a list of the value of the absorption (``\nu=+``) correlation function according to the given time list.
- `cmlist::Vector{ComplexF64}` : a list of the value of the emission (``\nu=-``) correlation function according to the given time list.
"""
function C(bath::FermionBath, tlist::AbstractVector)
    cplist = zeros(ComplexF64, length(tlist))
    cmlist = zeros(ComplexF64, length(tlist))
    for (i, t) in enumerate(tlist)
        for e in bath
            if e.types == "fA"
                cplist[i] += e.η * exp(- e.γ * t)
            else
                cmlist[i] += e.η * exp(- e.γ * t)
            end
        end
    end
    
    return cplist, cmlist
end
export BosonBath, BosonBathRWA
export AbstractBosonBath, bosonReal, bosonImag, bosonRealImag, bosonAbsorb, bosonEmit

abstract type AbstractBosonBath end

Base.show(io::IO, B::AbstractBosonBath) =
    print(io, "$(typeof(B))-type bath with $(B.Nterm) exponential-expansion terms\n")
Base.show(io::IO, m::MIME"text/plain", B::AbstractBosonBath) = show(io, B)

function _check_bosonic_coupling_operator(op)
    _op = HandleMatrixType(op, "op (coupling operator)", type = Operator)
    if !ishermitian(_op)
        @warn "The system-bosonic-bath coupling operator \"op\" should be Hermitian Operator."
    end
    return _op
end

_check_bosonic_RWA_coupling_operator(op) = HandleMatrixType(op, "op (coupling operator)", type = Operator)

function _combine_same_gamma(η::Vector{Ti}, γ::Vector{Tj}) where {Ti,Tj<:Number}
    if length(η) != length(γ)
        error("The length of \'η\' and \'γ\' should be the same.")
    end

    ηnew = ComplexF64[η[1]]
    γnew = ComplexF64[γ[1]]
    for j in 2:length(γ)
        valueDontExist = true
        for k in 1:length(γnew)
            if γnew[k] ≈ γ[j]
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
- `Nterm` : the number of exponential-expansion term of correlation functions
- `δ` : The approximation discrepancy which is used for adding the terminator to HEOM matrix (see function: addTerminator)

# Methods
One can obtain the ``k``-th exponent (exponential-expansion term) from `bath::BosonBath` by calling : `bath[k]`.
`HierarchicalEOM.jl` also supports the following calls (methods) :
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
    op::QuantumObject
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
- `op::QuantumObject` : The system coupling operator, must be Hermitian and, for fermionic systems, even-parity to be compatible with charge conservation.
- `η::Vector{Ti<:Number}` : the coefficients ``\eta_i`` in bath correlation function ``C(\tau)``.
- `γ::Vector{Tj<:Number}` : the coefficients ``\gamma_i`` in bath correlation function ``C(\tau)``.
- `δ::Number` : The approximation discrepancy (Default to `0.0`) which is used for adding the terminator to HEOM matrix (see function: addTerminator)
- `combine::Bool` : Whether to combine the exponential-expansion terms with the same frequency. Defaults to `true`.
"""
function BosonBath(
    op::QuantumObject,
    η::Vector{Ti},
    γ::Vector{Tj},
    δ::Number = 0.0;
    combine::Bool = true,
) where {Ti,Tj<:Number}
    _op = HandleMatrixType(op, "op (coupling operator)", type = Operator)
    if combine
        ηnew, γnew = _combine_same_gamma(η, γ)
        bRI = bosonRealImag(_op, real.(ηnew), imag.(ηnew), γnew)
    else
        bRI = bosonRealImag(_op, real.(η), imag.(η), γ)
    end
    return BosonBath(AbstractBosonBath[bRI], _op, bRI.Nterm, δ)
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
- `op::QuantumObject` : The system coupling operator, must be Hermitian and, for fermionic systems, even-parity to be compatible with charge conservation.
- `η_real::Vector{Ti<:Number}` : the coefficients ``\eta_i`` in real part of bath correlation function ``C^{u=\textrm{R}}``.
- `γ_real::Vector{Tj<:Number}` : the coefficients ``\gamma_i`` in real part of bath correlation function ``C^{u=\textrm{R}}``.
- `η_imag::Vector{Tk<:Number}` : the coefficients ``\eta_i`` in imaginary part of bath correlation function ``C^{u=\textrm{I}}``.
- `γ_imag::Vector{Tl<:Number}` : the coefficients ``\gamma_i`` in imaginary part of bath correlation function ``C^{u=\textrm{I}}``.
- `δ::Number` : The approximation discrepancy (Default to `0.0`) which is used for adding the terminator to HEOM matrix (see function: addTerminator)
- `combine::Bool` : Whether to combine the exponential-expansion terms with the same frequency. Defaults to `true`.
"""
function BosonBath(
    op::QuantumObject,
    η_real::Vector{Ti},
    γ_real::Vector{Tj},
    η_imag::Vector{Tk},
    γ_imag::Vector{Tl},
    δ::Tm = 0.0;
    combine::Bool = true,
) where {Ti,Tj,Tk,Tl,Tm<:Number}
    _op = HandleMatrixType(op, "op (coupling operator)", type = Operator)

    if combine
        ηR, γR = _combine_same_gamma(η_real, γ_real)
        ηI, γI = _combine_same_gamma(η_imag, γ_imag)

        NR = length(γR)
        NI = length(γI)
        Nterm = NR + NI

        ηRI_real = ComplexF64[]
        ηRI_imag = ComplexF64[]
        γRI = ComplexF64[]
        real_idx = Int64[]
        imag_idx = Int64[]
        for j in 1:NR
            for k in 1:NI
                if γR[j] ≈ γI[k]
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

        bRI = bosonRealImag(_op, ηRI_real, ηRI_imag, γRI)
        bR = bosonReal(_op, ηR, γR)
        bI = bosonImag(_op, ηI, γI)
        Nterm_new = bRI.Nterm + bR.Nterm + bI.Nterm
        if Nterm != (Nterm_new + bRI.Nterm)
            error("Conflicts occur in combining real and imaginary parts of bath correlation function.")
        end

        return BosonBath(AbstractBosonBath[bRI, bR, bI], _op, Nterm_new, δ)

    else
        bR = bosonReal(_op, η_real, γ_real)
        bI = bosonImag(_op, η_imag, γ_imag)
        return BosonBath(AbstractBosonBath[bR, bI], _op, bR.Nterm + bI.Nterm, δ)
    end
end

@doc raw"""
    struct bosonReal <: AbstractBosonBath
A bosonic bath for the real part of bath correlation function ``C^{u=\textrm{R}}``

# Fields
- `Comm`  : the super-operator (commutator) for the coupling operator.
- `dimensions` : the dimension list of the coupling operator (should be equal to the system dimensions).
- `η` : the coefficients ``\eta_i`` in real part of bath correlation function ``C^{u=\textrm{R}}``.
- `γ` : the coefficients ``\gamma_i`` in real part of bath correlation function ``C^{u=\textrm{R}}``.
- `Nterm` : the number of exponential-expansion term of correlation function

!!! note "`dims` property"
    For a given `b::bosonReal`, `b.dims` or `getproperty(b, :dims)` returns its `dimensions` in the type of integer-vector.
"""
struct bosonReal <: AbstractBosonBath
    Comm::SparseMatrixCSC{ComplexF64,Int64}
    dimensions::Dimensions
    η::AbstractVector
    γ::AbstractVector
    Nterm::Int
end

@doc raw"""
    bosonReal(op, η_real, γ_real)
Generate bosonic bath for the real part of bath correlation function ``C^{u=\textrm{R}}``

# Parameters
- `op::QuantumObject` : The system coupling operator, must be Hermitian and, for fermionic systems, even-parity to be compatible with charge conservation.
- `η_real::Vector{Ti<:Number}` : the coefficients ``\eta_i`` in real part of bath correlation function ``C^{u=\textrm{R}}``.
- `γ_real::Vector{Tj<:Number}` : the coefficients ``\gamma_i`` in real part of bath correlation function ``C^{u=\textrm{R}}``.
"""
function bosonReal(op::QuantumObject, η_real::Vector{Ti}, γ_real::Vector{Tj}) where {Ti,Tj<:Number}
    _op = _check_bosonic_coupling_operator(op)

    # check if the length of coefficients are valid
    N_exp_term = length(η_real)
    if N_exp_term != length(γ_real)
        error("The length of \'η_real\' and \'γ_real\' should be the same.")
    end
    Id_cache = I(size(_op, 1))
    return bosonReal(_spre(_op.data, Id_cache) - _spost(_op.data, Id_cache), _op.dimensions, η_real, γ_real, N_exp_term)
end

@doc raw"""
    struct bosonImag <: AbstractBosonBath
A bosonic bath for the imaginary part of bath correlation function ``C^{u=\textrm{I}}``

# Fields
- `Comm`  : the super-operator (commutator) for the coupling operator.
- `anComm`  : the super-operator (anti-commutator) for the coupling operator.
- `dimensions` : the dimension list of the coupling operator (should be equal to the system dimensions).
- `η` : the coefficients ``\eta_i`` in imaginary part of bath correlation function ``C^{u=\textrm{I}}``.
- `γ` : the coefficients ``\gamma_i`` in imaginary part of bath correlation function ``C^{u=\textrm{I}}``.
- `Nterm` : the number of exponential-expansion term of correlation function

!!! note "`dims` property"
    For a given `b::bosonImag`, `b.dims` or `getproperty(b, :dims)` returns its `dimensions` in the type of integer-vector.
"""
struct bosonImag <: AbstractBosonBath
    Comm::SparseMatrixCSC{ComplexF64,Int64}
    anComm::SparseMatrixCSC{ComplexF64,Int64}
    dimensions::Dimensions
    η::AbstractVector
    γ::AbstractVector
    Nterm::Int
end

@doc raw"""
    bosonImag(op, η_imag, γ_imag)
Generate bosonic bath for the imaginary part of correlation function ``C^{u=\textrm{I}}``

# Parameters
- `op::QuantumObject` : The system coupling operator, must be Hermitian and, for fermionic systems, even-parity to be compatible with charge conservation.
- `η_imag::Vector{Ti<:Number}` : the coefficients ``\eta_i`` in imaginary part of bath correlation functions ``C^{u=\textrm{I}}``.
- `γ_imag::Vector{Tj<:Number}` : the coefficients ``\gamma_i`` in imaginary part of bath correlation functions ``C^{u=\textrm{I}}``.
"""
function bosonImag(op::QuantumObject, η_imag::Vector{Ti}, γ_imag::Vector{Tj}) where {Ti,Tj<:Number}
    _op = _check_bosonic_coupling_operator(op)

    # check if the length of coefficients are valid
    N_exp_term = length(η_imag)
    if N_exp_term != length(γ_imag)
        error("The length of \'η_imag\' and \'γ_imag\' should be the same.")
    end
    Id_cache = I(size(_op, 1))
    spreQ = _spre(_op.data, Id_cache)
    spostQ = _spost(_op.data, Id_cache)
    return bosonImag(spreQ - spostQ, spreQ + spostQ, _op.dimensions, η_imag, γ_imag, N_exp_term)
end

@doc raw"""
    struct bosonRealImag <: AbstractBosonBath
A bosonic bath which the real part and imaginary part of the bath correlation function are combined 

# Fields
- `Comm`  : the super-operator (commutator) for the coupling operator.
- `anComm`  : the super-operator (anti-commutator) for the coupling operator.
- `dimensions` : the dimension list of the coupling operator (should be equal to the system dimensions).
- `η_real` : the real part of coefficients ``\eta_i`` in bath correlation function ``\sum_i \eta_i \exp(-\gamma_i t)``.
- `η_imag` : the imaginary part of coefficients ``\eta_i`` in bath correlation function ``\sum_i \eta_i \exp(-\gamma_i t)``.
- `γ` : the coefficients ``\gamma_i`` in bath correlation function ``\sum_i \eta_i \exp(-\gamma_i t)``.
- `Nterm` : the number of exponential-expansion term of correlation function

!!! note "`dims` property"
    For a given `b::bosonRealImag`, `b.dims` or `getproperty(b, :dims)` returns its `dimensions` in the type of integer-vector.
"""
struct bosonRealImag <: AbstractBosonBath
    Comm::SparseMatrixCSC{ComplexF64,Int64}
    anComm::SparseMatrixCSC{ComplexF64,Int64}
    dimensions::Dimensions
    η_real::AbstractVector
    η_imag::AbstractVector
    γ::AbstractVector
    Nterm::Int
end

@doc raw"""
    bosonRealImag(op, η_real, η_imag, γ)
Generate bosonic bath which the real part and imaginary part of the bath correlation function are combined

# Parameters
- `op::QuantumObject` : The system coupling operator, must be Hermitian and, for fermionic systems, even-parity to be compatible with charge conservation.
- `η_real::Vector{Ti<:Number}` : the real part of coefficients ``\eta_i`` in bath correlation function ``\sum_i \eta_i \exp(-\gamma_i t)``.
- `η_imag::Vector{Tj<:Number}` : the imaginary part of coefficients ``\eta_i`` in bath correlation function ``\sum_i \eta_i \exp(-\gamma_i t)``.
- `γ::Vector{Tk<:Number}` : the coefficients ``\gamma_i`` in bath correlation function ``\sum_i \eta_i \exp(-\gamma_i t)``.
"""
function bosonRealImag(
    op::QuantumObject,
    η_real::Vector{Ti},
    η_imag::Vector{Tj},
    γ::Vector{Tk},
) where {Ti,Tj,Tk<:Number}
    _op = _check_bosonic_coupling_operator(op)

    # check if the length of coefficients are valid
    N_exp_term = length(η_real)
    if (N_exp_term != length(η_imag)) || (N_exp_term != length(γ))
        error("The length of \'η_real\', \'η_imag\' and \'γ\' should be the same.")
    end
    Id_cache = I(size(_op, 1))
    spreQ = _spre(_op.data, Id_cache)
    spostQ = _spost(_op.data, Id_cache)
    return bosonRealImag(spreQ - spostQ, spreQ + spostQ, _op.dimensions, η_real, η_imag, γ, N_exp_term)
end

@doc raw"""
    BosonBathRWA(op, η_absorb, γ_absorb, η_emit, γ_emit, δ=0.0)
A function for generating BosonBath object where the interaction between system and bosonic bath applies the rotating wave approximation (RWA).

```math
\begin{aligned}
C^{\nu=+}(\tau)
&=\frac{1}{2\pi}\int_{0}^{\infty} d\omega J(\omega) n(\omega) e^{i\omega \tau}\\
&=\sum_i \eta_i^{\nu=+} \exp(-\gamma_i^{\nu=+} \tau),\\
C^{\nu=-}(\tau)
&=\frac{1}{2\pi}\int_{0}^{\infty} d\omega J(\omega) (1+n(\omega)) e^{-i\omega \tau}\\
&=\sum_i \eta_i^{\nu=-} \exp(-\gamma_i^{\nu=-} \tau),
\end{aligned}
```
where ``\nu=+`` (``\nu=-``) represents absorption (emission) process, ``J(\omega)`` is the spectral density of the bath and ``n(\omega)`` is the Bose-Einstein distribution.

# Parameters
- `op::QuantumObject` : The system annihilation operator according to the system-bosonic-bath interaction.
- `η_absorb::Vector{Ti<:Number}` : the coefficients ``\eta_i`` of absorption bath correlation function ``C^{\nu=+}(\tau)``.
- `γ_absorb::Vector{Tj<:Number}` : the coefficients ``\gamma_i`` of absorption bath correlation function ``C^{\nu=+}(\tau)``.
- `η_emit::Vector{Tk<:Number}` : the coefficients ``\eta_i`` of emission bath correlation function ``C^{\nu=-}(\tau)``.
- `γ_emit::Vector{Tl<:Number}` : the coefficients ``\gamma_i`` of emission bath correlation function ``C^{\nu=-}(\tau)``.
- `δ::Number` : The approximation discrepancy (Defaults to `0.0`) which is used for adding the terminator to HEOMLS matrix (see function: addTerminator)
"""
function BosonBathRWA(
    op::QuantumObject,
    η_absorb::Vector{Ti},
    γ_absorb::Vector{Tj},
    η_emit::Vector{Tk},
    γ_emit::Vector{Tl},
    δ::Tm = 0.0,
) where {Ti,Tj,Tk,Tl,Tm<:Number}
    _check_gamma_absorb_and_emit(γ_absorb, γ_emit)
    _op = HandleMatrixType(op, "op (coupling operator)", type = Operator)

    bA = bosonAbsorb(adjoint(_op), η_absorb, γ_absorb, η_emit)
    bE = bosonEmit(_op, η_emit, γ_emit, η_absorb)
    return BosonBath(AbstractBosonBath[bA, bE], _op, bA.Nterm + bE.Nterm, δ)
end

@doc raw"""
    struct bosonAbsorb <: AbstractBosonBath
An bath object which describes the absorption process of the bosonic system by a correlation function ``C^{\nu=+}``

# Fields
- `spre`   : the super-operator (left side operator multiplication) for the coupling operator.
- `spost`  : the super-operator (right side operator multiplication) for the coupling operator.
- `CommD`  : the super-operator (commutator) for the adjoint of the coupling operator.
- `dimensions` : the dimension list of the coupling operator (should be equal to the system dimensions).
- `η` : the coefficients ``\eta_i`` of absorption bath correlation function ``C^{\nu=+}``.
- `γ` : the coefficients ``\gamma_i`` of absorption bath correlation function ``C^{\nu=+}``.
- `η_emit` : the coefficients ``\eta_i`` of emission bath correlation function ``C^{\nu=-}``.
- `Nterm` : the number of exponential-expansion term of correlation function

!!! note "`dims` property"
    For a given `b::bosonAbsorb`, `b.dims` or `getproperty(b, :dims)` returns its `dimensions` in the type of integer-vector.
"""
struct bosonAbsorb <: AbstractBosonBath
    spre::SparseMatrixCSC{ComplexF64,Int64}
    spost::SparseMatrixCSC{ComplexF64,Int64}
    CommD::SparseMatrixCSC{ComplexF64,Int64}
    dimensions::Dimensions
    η::AbstractVector
    γ::AbstractVector
    η_emit::AbstractVector
    Nterm::Int
end

@doc raw"""
    bosonAbsorb(op, η_absorb, γ_absorb, η_emit)
Generate bosonic bath which describes the absorption process of the bosonic system by a correlation function ``C^{\nu=+}``

# Parameters
- `op::QuantumObject` : The system creation operator according to the system-fermionic-bath interaction.
- `η_absorb::Vector{Ti<:Number}` : the coefficients ``\eta_i`` of absorption bath correlation function ``C^{\nu=+}``.
- `γ_absorb::Vector{Tj<:Number}` : the coefficients ``\gamma_i`` of absorption bath correlation function ``C^{\nu=+}``.
- `η_emit::Vector{Tk<:Number}` : the coefficients ``\eta_i`` of emission bath correlation function ``C^{\nu=-}``.
"""
function bosonAbsorb(
    op::QuantumObject,
    η_absorb::Vector{Ti},
    γ_absorb::Vector{Tj},
    η_emit::Vector{Tk},
) where {Ti,Tj,Tk<:Number}
    _op = _check_bosonic_RWA_coupling_operator(op)

    # check if the length of coefficients are valid
    N_exp_term = length(η_absorb)
    if (N_exp_term != length(γ_absorb)) || (N_exp_term != length(η_emit))
        error("The length of \'η_absorb\', \'γ_absorb\' and \'η_emit\' should all be the same.")
    end
    Id_cache = I(size(_op, 1))
    return bosonAbsorb(
        _spre(_op.data, Id_cache),
        _spost(_op.data, Id_cache),
        _spre(adjoint(_op).data, Id_cache) - _spost(adjoint(_op).data, Id_cache),
        _op.dimensions,
        η_absorb,
        γ_absorb,
        η_emit,
        N_exp_term,
    )
end

@doc raw"""
    struct bosonEmit <: AbstractBosonBath
An bath object which describes the emission process of the bosonic system by a correlation function ``C^{\nu=-}``

# Fields
- `spre`   : the super-operator (left side operator multiplication) for the coupling operator.
- `spost`  : the super-operator (right side operator multiplication) for the coupling operator.
- `CommD`  : the super-operator (commutator) for the adjoint of the coupling operator.
- `dimensions` : the dimension list of the coupling operator (should be equal to the system dimensions).
- `η` : the coefficients ``\eta_i`` of emission bath correlation function ``C^{\nu=-}``.
- `γ` : the coefficients ``\gamma_i`` of emission bath correlation function ``C^{\nu=-}``.
- `η_absorb` : the coefficients ``\eta_i`` of absorption bath correlation function ``C^{\nu=+}``.
- `Nterm` : the number of exponential-expansion term of correlation function

!!! note "`dims` property"
    For a given `b::bosonEmit`, `b.dims` or `getproperty(b, :dims)` returns its `dimensions` in the type of integer-vector.
"""
struct bosonEmit <: AbstractBosonBath
    spre::SparseMatrixCSC{ComplexF64,Int64}
    spost::SparseMatrixCSC{ComplexF64,Int64}
    CommD::SparseMatrixCSC{ComplexF64,Int64}
    dimensions::Dimensions
    η::AbstractVector
    γ::AbstractVector
    η_absorb::AbstractVector
    Nterm::Int
end

@doc raw"""
    bosonEmit(op, η_emit, γ_emit, η_absorb)
Generate bosonic bath which describes the emission process of the bosonic system by a correlation function ``C^{\nu=-}``

# Parameters
- `op::QuantumObject` : The system annihilation operator according to the system-bosonic-bath interaction.
- `η_emit::Vector{Ti<:Number}` : the coefficients ``\eta_i`` of emission bath correlation function ``C^{\nu=-}``.
- `γ_emit::Vector{Ti<:Number}` : the coefficients ``\gamma_i`` of emission bath correlation function ``C^{\nu=-}``.
- `η_absorb::Vector{Ti<:Number}` : the coefficients ``\eta_i`` of absorption bath correlation function ``C^{\nu=+}``.
"""
function bosonEmit(
    op::QuantumObject,
    η_emit::Vector{Ti},
    γ_emit::Vector{Tj},
    η_absorb::Vector{Tk},
) where {Ti,Tj,Tk<:Number}
    _op = _check_bosonic_RWA_coupling_operator(op)

    # check if the length of coefficients are valid
    N_exp_term = length(η_emit)
    if (N_exp_term != length(γ_emit)) || (N_exp_term != length(η_absorb))
        error("The length of \'η_emit\', \'γ_emit\' and \'η_absorb\' should all be the same.")
    end
    Id_cache = I(size(_op, 1))
    return bosonEmit(
        _spre(_op.data, Id_cache),
        _spost(_op.data, Id_cache),
        _spre(adjoint(_op).data, Id_cache) - _spost(adjoint(_op).data, Id_cache),
        _op.dimensions,
        η_emit,
        γ_emit,
        η_absorb,
        N_exp_term,
    )
end

function Base.getproperty(
    b::BType,
    key::Symbol,
) where {B<:Union{bosonReal,bosonImage,bosonRealImag,bosonAbsorb,bosonEmit}}
    # a comment here to avoid bad render by JuliaFormatter
    if key === :dims
        return dimensions_to_dims(getfield(b, :dimensions))
    else
        return getfield(b, key)
    end
end

@doc raw"""
    correlation_function(bath, tlist)
Calculate the correlation function ``C(t)`` for a given bosonic bath and time list.

## if the input bosonic bath did not apply rotating wave approximation (RWA)
```math
C(t)=\sum_{u=\textrm{R},\textrm{I}}(\delta_{u, \textrm{R}} + i\delta_{u, \textrm{I}})C^{u}(t)
```
where
```math
C^{u}(t)=\sum_i \eta_i^u e^{-\gamma_i^u t}
```

## if the input bosonic bath applies rotating wave approximation (RWA)
```math
C^{\nu=\pm}(t)=\sum_i \eta_i^\nu e^{-\gamma_i^\nu t}
```

# Parameters
- `bath::BosonBath` : The bath object which describes a certain bosonic bath.
- `tlist::AbstractVector`: The specific time.

# Returns (without RWA)
- `clist::Vector{ComplexF64}` : a list of the value of correlation function according to the given time list.

# Returns (with RWA)
- `cplist::Vector{ComplexF64}` : a list of the value of the absorption (``\nu=+``) correlation function according to the given time list.
- `cmlist::Vector{ComplexF64}` : a list of the value of the emission (``\nu=-``) correlation function according to the given time list.
"""
function correlation_function(bath::BosonBath, tlist::AbstractVector)
    T = (bath[1]).types

    # without RWA
    if (T == "bR") || (T == "bI") || (T == "bRI")
        clist = zeros(ComplexF64, length(tlist))
        for (i, t) in enumerate(tlist)
            for e in bath
                if (e.types == "bR") || (e.types == "bRI")
                    clist[i] += e.η * exp(-e.γ * t)
                elseif e.types == "bI"
                    clist[i] += 1.0im * e.η * exp(-e.γ * t)
                else
                    error("Invalid bath")
                end
            end
        end
        return clist

    else # with RWA
        cplist = zeros(ComplexF64, length(tlist))
        cmlist = zeros(ComplexF64, length(tlist))
        for (i, t) in enumerate(tlist)
            for e in bath
                if e.types == "bA"
                    cplist[i] += e.η * exp(-e.γ * t)
                elseif e.types == "bE"
                    cmlist[i] += e.η * exp(-e.γ * t)
                else
                    error("Invalid bath")
                end
            end
        end
        return cplist, cmlist
    end
end

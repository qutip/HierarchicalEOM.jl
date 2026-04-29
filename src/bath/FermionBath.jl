export FermionBath
export AbstractFermionBath, fermionAbsorb, fermionEmit

abstract type AbstractFermionBath end

Base.show(io::IO, B::AbstractFermionBath) =
    print(io, "$(typeof(B))-type bath with $(B.Nterm) exponential-expansion terms\n")
Base.show(io::IO, m::MIME"text/plain", B::AbstractFermionBath) = show(io, B)

function _check_fermionic_coupling_operator(op)
    _op = HandleMatrixType(op, "op (coupling operator)", type = Operator())
    issparse(_op.data) || (@warn "The system-fermionic-bath coupling operator \"op\" is recommended to be a sparse matrix for better performance.")
    return _op
end

@doc raw"""
    struct FermionBath <: AbstractBath
An object which describes the interaction between system and fermionic bath

# Fields
- `bath` : the different fermion-bath-type objects which describes the interaction
- `op` : The system \"emission\" operator according to the system-fermionic-bath interaction.
- `Nterm` : the total number of different bath terms.
- `δ` : The approximation discrepancy which is used for adding the terminator to HEOM matrix (see function: addTerminator)

# Methods
One can obtain the ``k``-th term from `bath::FermionBath` by calling : `bath[k]`.
`HierarchicalEOM.jl` also supports the following calls (methods) :
```julia
bath[1:k];   # returns a vector which contains the `1`-st to the `k`-th term.
bath[1:end]; # returns a vector which contains all terms
bath[:];     # returns a vector which contains all terms
from b in bath
    # do something
end
```
"""
struct FermionBath <: AbstractBath
    bath::Vector{AbstractFermionBath}
    op::QuantumObject
    Nterm::Int
    δ::Number

    function FermionBath(bath::Vector{AbstractFermionBath}, op::QuantumObject, Nterm::Int, δ::Number)
        (Nterm == sum(getfield.(bath, :Nterm))) || error("Invalid `Nterm` in FermionBath.")
        return new(bath, op, Nterm, δ)
    end
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
- `op::QuantumObject` : The system annihilation operator according to the system-fermionic-bath interaction.
- `η_absorb::Vector{Ti<:Number}` : the coefficients ``\eta_i`` of absorption bath correlation function ``C^{\nu=+}(\tau)``.
- `γ_absorb::Vector{Tj<:Number}` : the coefficients ``\gamma_i`` of absorption bath correlation function ``C^{\nu=+}(\tau)``.
- `η_emit::Vector{Tk<:Number}` : the coefficients ``\eta_i`` of emission bath correlation function ``C^{\nu=-}(\tau)``.
- `γ_emit::Vector{Tl<:Number}` : the coefficients ``\gamma_i`` of emission bath correlation function ``C^{\nu=-}(\tau)``.
- `δ::Number` : The approximation discrepancy (Defaults to `0.0`) which is used for adding the terminator to HEOMLS matrix (see function: addTerminator)
"""
function FermionBath(
        op::QuantumObject,
        η_absorb::Vector{Ti},
        γ_absorb::Vector{Tj},
        η_emit::Vector{Tk},
        γ_emit::Vector{Tl},
        δ::Tm = 0.0,
    ) where {Ti <: Number, Tj <: Number, Tk <: Number, Tl <: Number, Tm <: Number}
    _check_gamma_absorb_and_emit(γ_absorb, γ_emit)

    _op = HandleMatrixType(op, "op (coupling operator)", type = Operator())
    fA = fermionAbsorb(adjoint(_op), η_absorb, γ_absorb, η_emit)
    fE = fermionEmit(_op, η_emit, γ_emit, η_absorb)
    return FermionBath(AbstractFermionBath[fA, fE], _op, fA.Nterm + fE.Nterm, δ)
end

@doc raw"""
    struct fermionAbsorb <: AbstractFermionBath
An bath object which describes the absorption process of the fermionic system by a correlation function ``C^{\nu=+}``

# Fields
- `spre`   : the super-operator (left side operator multiplication) for the coupling operator.
- `spost`  : the super-operator (right side operator multiplication) for the coupling operator.
- `spreD`  : the super-operator (left side operator multiplication) for the adjoint of the coupling operator.
- `spostD` : the super-operator (right side operator multiplication) for the adjoint of the coupling operator.
- `dimensions` : the dimension list of the coupling operator (should be equal to the system dimensions).
- `η` : the coefficients ``\eta_i`` of absorption bath correlation function ``C^{\nu=+}``.
- `γ` : the coefficients ``\gamma_i`` of absorption bath correlation function ``C^{\nu=+}``.
- `η_emit` : the coefficients ``\eta_i`` of emission bath correlation function ``C^{\nu=-}``.
- `Nterm` : the number of exponential-expansion term of correlation function.

!!! note "`dims` property"
    For a given `b::fermionAbsorb`, `b.dims` or `getproperty(b, :dims)` returns its `dimensions` in the type of integer-vector.
"""
struct fermionAbsorb <: AbstractFermionBath
    spre::SparseMatrixCSC{ComplexF64, Int64}
    spost::SparseMatrixCSC{ComplexF64, Int64}
    spreD::SparseMatrixCSC{ComplexF64, Int64}
    spostD::SparseMatrixCSC{ComplexF64, Int64}
    dimensions::Dimensions
    η::AbstractVector
    γ::AbstractVector
    η_emit::AbstractVector
    Nterm::Int
end

@doc raw"""
    fermionAbsorb(op, η_absorb, γ_absorb, η_emit)
Generate fermionic bath which describes the absorption process of the fermionic system by a correlation function ``C^{\nu=+}``

# Parameters
- `op::QuantumObject` : The system creation operator according to the system-fermionic-bath interaction.
- `η_absorb::Vector{Ti<:Number}` : the coefficients ``\eta_i`` of absorption bath correlation function ``C^{\nu=+}``.
- `γ_absorb::Vector{Tj<:Number}` : the coefficients ``\gamma_i`` of absorption bath correlation function ``C^{\nu=+}``.
- `η_emit::Vector{Tk<:Number}` : the coefficients ``\eta_i`` of emission bath correlation function ``C^{\nu=-}``.
"""
function fermionAbsorb(
        op::QuantumObject,
        η_absorb::Vector{Ti},
        γ_absorb::Vector{Tj},
        η_emit::Vector{Tk},
    ) where {Ti <: Number, Tj <: Number, Tk <: Number}
    _op = _check_fermionic_coupling_operator(op)

    # check if the length of coefficients are valid
    N_exp_term = length(η_absorb)
    if (N_exp_term != length(γ_absorb)) || (N_exp_term != length(η_emit))
        error("The length of \'η_absorb\', \'γ_absorb\' and \'η_emit\' should all be the same.")
    end
    return fermionAbsorb(
        _spre(_op.data),
        _spost(_op.data),
        _spre(adjoint(_op).data),
        _spost(adjoint(_op).data),
        _op.dimensions,
        η_absorb,
        γ_absorb,
        η_emit,
        N_exp_term,
    )
end

@doc raw"""
    struct fermionEmit <: AbstractFermionBath
An bath object which describes the emission process of the fermionic system by a correlation function ``C^{\nu=-}``

# Fields
- `spre`   : the super-operator (left side operator multiplication) for the coupling operator.
- `spost`  : the super-operator (right side operator multiplication) for the coupling operator.
- `spreD`  : the super-operator (left side operator multiplication) for the adjoint of the coupling operator.
- `spostD` : the super-operator (right side operator multiplication) for the adjoint of the coupling operator.
- `dimensions` : the dimension list of the coupling operator (should be equal to the system dimensions).
- `η` : the coefficients ``\eta_i`` of emission bath correlation function ``C^{\nu=-}``.
- `γ` : the coefficients ``\gamma_i`` of emission bath correlation function ``C^{\nu=-}``.
- `η_absorb` : the coefficients ``\eta_i`` of absorption bath correlation function ``C^{\nu=+}``.
- `Nterm` : the number of exponential-expansion term of correlation function.

!!! note "`dims` property"
    For a given `b::fermionEmit`, `b.dims` or `getproperty(b, :dims)` returns its `dimensions` in the type of integer-vector.
"""
struct fermionEmit <: AbstractFermionBath
    spre::SparseMatrixCSC{ComplexF64, Int64}
    spost::SparseMatrixCSC{ComplexF64, Int64}
    spreD::SparseMatrixCSC{ComplexF64, Int64}
    spostD::SparseMatrixCSC{ComplexF64, Int64}
    dimensions::Dimensions
    η::AbstractVector
    γ::AbstractVector
    η_absorb::AbstractVector
    Nterm::Int
end

@doc raw"""
    fermionEmit(op, η_emit, γ_emit, η_absorb)
Generate fermionic bath which describes the emission process of the fermionic system by a correlation function ``C^{\nu=-}``

# Parameters
- `op::QuantumObject` : The system annihilation operator according to the system-fermionic-bath interaction.
- `η_emit::Vector{Ti<:Number}` : the coefficients ``\eta_i`` of emission bath correlation function ``C^{\nu=-}``.
- `γ_emit::Vector{Ti<:Number}` : the coefficients ``\gamma_i`` of emission bath correlation function ``C^{\nu=-}``.
- `η_absorb::Vector{Ti<:Number}` : the coefficients ``\eta_i`` of absorption bath correlation function ``C^{\nu=+}``.
"""
function fermionEmit(
        op::QuantumObject,
        η_emit::Vector{Ti},
        γ_emit::Vector{Tj},
        η_absorb::Vector{Tk},
    ) where {Ti <: Number, Tj <: Number, Tk <: Number}
    _op = _check_fermionic_coupling_operator(op)

    # check if the length of coefficients are valid
    N_exp_term = length(η_emit)
    if (N_exp_term != length(γ_emit)) || (N_exp_term != length(η_absorb))
        error("The length of \'η_emit\', \'γ_emit\' and \'η_absorb\' should all be the same.")
    end

    return fermionEmit(
        _spre(_op.data),
        _spost(_op.data),
        _spre(adjoint(_op).data),
        _spost(adjoint(_op).data),
        _op.dimensions,
        η_emit,
        γ_emit,
        η_absorb,
        N_exp_term,
    )
end

function Base.getproperty(b::BType, key::Symbol) where {BType <: Union{fermionAbsorb, fermionEmit}}
    # a comment here to avoid bad render by JuliaFormatter
    if key === :dims
        return dimensions_to_dims(getfield(b, :dimensions))
    else
        return getfield(b, key)
    end
end

@doc raw"""
    correlation_function(bath, tlist)
Calculate the correlation function ``C^{\nu=+}(t)`` and ``C^{\nu=-}(t)`` for a given fermionic bath and time list.
Here, ``\nu=+`` represents the absorption process and ``\nu=-`` represents the emission process.

```math
C^{\nu=\pm}(t)=\sum_i \eta_i^\nu e^{-\gamma_i^\nu t}
```

# Parameters
- `bath::FermionBath` : The bath object which describes a certain fermionic bath.
- `tlist::AbstractVector`: The specific time.

# Returns
- `cplist::Vector{ComplexF64}` : a list of the value of the absorption (``\nu=+``) correlation function according to the given time list.
- `cmlist::Vector{ComplexF64}` : a list of the value of the emission (``\nu=-``) correlation function according to the given time list.
"""
function correlation_function(bath::FermionBath, tlist::AbstractVector)
    cplist = zeros(ComplexF64, length(tlist))
    cmlist = zeros(ComplexF64, length(tlist))
    for (i, t) in enumerate(tlist)
        for e in bath
            if e.types == "fA"
                cplist[i] += e.η * exp(-e.γ * t)
            else
                cmlist[i] += e.η * exp(-e.γ * t)
            end
        end
    end

    return cplist, cmlist
end

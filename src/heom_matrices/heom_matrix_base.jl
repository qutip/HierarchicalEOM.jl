export HEOMSuperOp
export Propagator
export addBosonDissipator, addFermionDissipator, addTerminator
export get_cached_HEOMLS_data

@doc raw"""
    struct HEOMSuperOp
General HEOM superoperator matrix.  

# Fields
- `data` : the HEOM superoperator matrix
- `dimensions` : the dimension list of the coupling operator (should be equal to the system dimensions).
- `N` : the number of auxiliary density operators
- `parity`: the parity label (`EVEN` or `ODD`).

!!! note "`dims` property"
    For a given `M::HEOMSuperOp`, `M.dims` or `getproperty(M, :dims)` returns its `dimensions` in the type of integer-vector.
"""
struct HEOMSuperOp
    data::SparseMatrixCSC{ComplexF64,Int64}
    dimensions::Dimensions
    N::Int
    parity::AbstractParity
end

function Base.getproperty(M::HEOMSuperOp, key::Symbol)
    # a comment here to avoid bad render by JuliaFormatter
    if key === :dims
        return dimensions_to_dims(getfield(M, :dimensions))
    else
        return getfield(M, key)
    end
end

@doc raw"""
    HEOMSuperOp(op, opParity, refHEOMLS)
Construct the HEOM superoperator matrix corresponding to the given system SuperOperator which acts on all `ADOs`.  

During the multiplication on all the `ADOs`, the parity of the output `ADOs` might change depend on the parity of this HEOM superoperator.

# Parameters
- `op` : The system SuperOperator which will act on all `ADOs`.
- `opParity::AbstractParity` : the parity label of the given operator (`op`), should be `EVEN` or `ODD`.
- `refHEOMLS::AbstractHEOMLSMatrix` : copy the system `dimensions` and number of `ADOs` (`N`) from this reference HEOMLS matrix
"""
HEOMSuperOp(op, opParity::AbstractParity, refHEOMLS::AbstractHEOMLSMatrix) =
    HEOMSuperOp(op, opParity, refHEOMLS.dimensions, refHEOMLS.N)

@doc raw"""
    HEOMSuperOp(op, opParity, refADOs)
Construct the HEOM superoperator matrix corresponding to the given system SuperOperator which acts on all `ADOs`.  

During the multiplication on all the `ADOs`, the parity of the output `ADOs` might change depend on the parity of this HEOM superoperator.

# Parameters
- `op` : The system SuperOperator which will act on all `ADOs`.
- `opParity::AbstractParity` : the parity label of the given operator (`op`), should be `EVEN` or `ODD`.
- `refADOs::ADOs` : copy the system `dimensions` and number of `ADOs` (`N`) from this reference `ADOs`   
"""
HEOMSuperOp(op, opParity::AbstractParity, refADOs::ADOs) = HEOMSuperOp(op, opParity, refADOs.dimensions, refADOs.N)

@doc raw"""
    HEOMSuperOp(op, opParity, dims, N)
Construct the HEOM superoperator matrix corresponding to the given system SuperOperator which acts on all `ADOs`.  

During the multiplication on all the `ADOs`, the parity of the output `ADOs` might change depend on the parity of this HEOM superoperator.

# Parameters
- `op` : The system SuperOperator which will act on all `ADOs`.
- `opParity::AbstractParity` : the parity label of the given operator (`op`), should be `EVEN` or `ODD`.
- `dims` : the dimension list of the coupling operator (should be equal to the system `dimensions`).
- `N::Int` : the number of `ADOs`.
"""
function HEOMSuperOp(op, opParity::AbstractParity, dims, N::Int)
    dimensions = _gen_dimensions(dims)
    sup_op = HandleMatrixType(op, dimensions, "op (operator)"; type = SuperOperator())

    return HEOMSuperOp(kron(Eye(N), sup_op.data), dimensions, N, opParity)
end

@doc raw"""
    size(M::HEOMSuperOp)
Returns the size of the HEOM superoperator matrix
"""
Base.size(M::HEOMSuperOp) = size(M.data)

@doc raw"""
    size(M::HEOMSuperOp, dim::Int)
Returns the specified dimension of the HEOM superoperator matrix
"""
Base.size(M::HEOMSuperOp, dim::Int) = size(M.data, dim)

@doc raw"""
    size(M::AbstractHEOMLSMatrix)
Returns the size of the HEOM Liouvillian superoperator matrix
"""
Base.size(M::AbstractHEOMLSMatrix) = size(M.data)

@doc raw"""
    size(M::AbstractHEOMLSMatrix, dim::Int)
Returns the specified dimension of the HEOM Liouvillian superoperator matrix
"""
Base.size(M::AbstractHEOMLSMatrix, dim::Int) = size(M.data, dim)

@doc raw"""
    eltype(M::HEOMSuperOp)
Returns the elements' type of the HEOM superoperator matrix
"""
Base.eltype(M::HEOMSuperOp) = eltype(M.data)

@doc raw"""
    eltype(M::AbstractHEOMLSMatrix)
Returns the elements' type of the HEOM Liouvillian superoperator matrix
"""
Base.eltype(M::AbstractHEOMLSMatrix) = eltype(M.data)

Base.getindex(M::HEOMSuperOp, i::Ti, j::Tj) where {Ti,Tj<:Any} = M.data[i, j]
Base.getindex(M::AbstractHEOMLSMatrix, i::Ti, j::Tj) where {Ti,Tj<:Any} = M.data[i, j]

function Base.show(io::IO, M::HEOMSuperOp)
    print(
        io,
        "$(M.parity) HEOM superoperator matrix acting on arbitrary-parity-ADOs\n",
        "system dims = $(_get_dims_string(M.dimensions))\n",
        "number of ADOs N = $(M.N)\n",
        "data =\n",
    )
    return show(io, MIME("text/plain"), M.data)
end

function Base.show(io::IO, M::AbstractHEOMLSMatrix)
    T = typeof(M)
    if T <: M_S
        type = "Schrodinger Eq."
    elseif T <: M_Boson
        type = "Boson"
    elseif T <: M_Fermion
        type = "Fermion"
    else
        type = "Boson-Fermion"
    end

    print(
        io,
        type,
        " type HEOMLS matrix acting on $(M.parity) ADOs\n",
        "system dims = $(_get_dims_string(M.dimensions))\n",
        "number of ADOs N = $(M.N)\n",
        "data =\n",
    )
    return show(io, MIME("text/plain"), M.data)
end

Base.show(io::IO, m::MIME"text/plain", M::HEOMSuperOp) = show(io, M)
Base.show(io::IO, m::MIME"text/plain", M::AbstractHEOMLSMatrix) = show(io, M)

function Base.:(*)(Sup::HEOMSuperOp, ados::ADOs)
    _check_sys_dim_and_ADOs_num(Sup, ados)

    return ADOs(Sup.data * ados.data, ados.dimensions, ados.N, Sup.parity * ados.parity)
end

function Base.:(*)(Sup1::HEOMSuperOp, Sup2::HEOMSuperOp)
    _check_sys_dim_and_ADOs_num(Sup1, Sup2)

    return HEOMSuperOp(Sup1.data * Sup2.data, Sup1.dimensions, Sup1.N, Sup1.parity * Sup2.parity)
end

Base.:(*)(n::Number, Sup::HEOMSuperOp) = HEOMSuperOp(n * Sup.data, Sup.dimensions, Sup.N, Sup.parity)
Base.:(*)(Sup::HEOMSuperOp, n::Number) = n * Sup

function Base.:(+)(Sup1::HEOMSuperOp, Sup2::HEOMSuperOp)
    _check_sys_dim_and_ADOs_num(Sup1, Sup2)
    _check_parity(Sup1, Sup2)

    return HEOMSuperOp(Sup1.data + Sup2.data, Sup1.dimensions, Sup1.N, Sup1.parity)
end

Base.:(+)(M::MatrixOperator, Sup::HEOMSuperOp) = MatrixOperator(M.A + Sup.data)
Base.:(+)(M::AbstractSciMLOperator, Sup::HEOMSuperOp) = M + Sup.data

function Base.:(+)(M::AbstractHEOMLSMatrix, Sup::HEOMSuperOp)
    _check_sys_dim_and_ADOs_num(M, Sup)
    _check_parity(M, Sup)

    return _reset_HEOMLS_data(M, M.data + Sup)
end

function Base.:(-)(Sup1::HEOMSuperOp, Sup2::HEOMSuperOp)
    _check_sys_dim_and_ADOs_num(Sup1, Sup2)
    _check_parity(Sup1, Sup2)

    return HEOMSuperOp(Sup1.data - Sup2.data, Sup1.dimensions, Sup1.N, Sup1.parity)
end

Base.:(-)(M::MatrixOperator, Sup::HEOMSuperOp) = MatrixOperator(M.A - Sup.data)
Base.:(-)(M::AbstractSciMLOperator, Sup::HEOMSuperOp) = M - Sup.data

function Base.:(-)(M::AbstractHEOMLSMatrix, Sup::HEOMSuperOp)
    _check_sys_dim_and_ADOs_num(M, Sup)
    _check_parity(M, Sup)

    return _reset_HEOMLS_data(M, M.data - Sup)
end

get_cached_HEOMLS_data(M::AbstractHEOMLSMatrix, cachevec::AbstractVector) = get_cached_HEOMLS_data(M.data, cachevec)

function get_cached_HEOMLS_data(M::T, cachevec::AbstractVector) where {T<:SciMLOperators.AddedOperator}
    ops = M.ops
    tensor_cache = nothing

    cached_op = sum(ops) do op
        if op isa TensorProductOperator
            tensor_cache === nothing && (tensor_cache = SciMLOperators.cache_operator(op, cachevec).cache)
            TensorProductOperator(op.ops, tensor_cache)
        else
            cache_operator(op, cachevec)
        end
    end

    return cached_op
end

get_cached_HEOMLS_data(M::T, cachevec::AbstractVector) where {T<:SciMLOperators.MatrixOperator} = M

@doc raw"""
    SciMLOperators.iscached(M::AbstractHEOMLSMatrix)

Test whether the [`AbstractHEOMLSMatrix`](@ref) `M` has preallocated caches for inplace evaluations.
"""
SciMLOperators.iscached(M::AbstractHEOMLSMatrix) = iscached(M.data)

@doc raw"""
    SciMLOperators.isconstant(M::AbstractHEOMLSMatrix)

Test whether the [`AbstractHEOMLSMatrix`](@ref) `M` is constant in time.
"""
SciMLOperators.isconstant(M::AbstractHEOMLSMatrix) = isconstant(M.data)

@doc raw"""
    SciMLOperators.concretize(M::AbstractHEOMLSMatrix)

Convert `M` to a concrete matrix representation by evaluating all lazy operations.
"""
SciMLOperators.concretize(M::AbstractHEOMLSMatrix) = _reset_HEOMLS_data(M, MatrixOperator(concretize(M.data)))

@doc raw"""
    Propagator(M, Δt; threshold, nonzero_tol)
Use `FastExpm.jl` to calculate the propagator matrix from a given HEOM Liouvillian superoperator matrix ``M`` with a specific time step ``\Delta t``.
That is, ``\exp(M * \Delta t)``.

# Parameters
- `M::AbstractHEOMLSMatrix` : the matrix given from HEOM model
- `Δt::Real` : A specific time step (time interval).
- `threshold::Real` : Determines the threshold for the Taylor series. Defaults to `1.0e-6`.
- `nonzero_tol::Real` : Strips elements smaller than `nonzero_tol` at each computation step to preserve sparsity. Defaults to `1.0e-14`.

For more details, please refer to [`FastExpm.jl`](https://github.com/fmentink/FastExpm.jl)

# Returns
- `::SparseMatrixCSC{ComplexF64, Int64}` : the propagator matrix
"""
@noinline function Propagator(
    M::AbstractHEOMLSMatrix{<:MatrixOperator},
    Δt::Real;
    threshold = 1.0e-6,
    nonzero_tol = 1.0e-14,
)
    return fastExpm(M.data.A * Δt; threshold = threshold, nonzero_tol = nonzero_tol)
end

function _reset_HEOMLS_data(M::T, new_data::AbstractSciMLOperator) where {T<:AbstractHEOMLSMatrix}
    if T <: M_S
        return M_S(new_data, M.tier, M.dimensions, M.N, M.sup_dim, M.parity)
    elseif T <: M_Boson
        return M_Boson(new_data, M.tier, M.dimensions, M.N, M.sup_dim, M.parity, M.bath, M.hierarchy)
    elseif T <: M_Fermion
        return M_Fermion(new_data, M.tier, M.dimensions, M.N, M.sup_dim, M.parity, M.bath, M.hierarchy)
    else
        return M_Boson_Fermion(
            new_data,
            M.Btier,
            M.Ftier,
            M.dimensions,
            M.N,
            M.sup_dim,
            M.parity,
            M.Bbath,
            M.Fbath,
            M.hierarchy,
        )
    end
end

@doc raw"""
    addBosonDissipator(M, jumpOP)
Adding bosonic dissipator to a given HEOMLS matrix which describes how the system dissipatively interacts with an extra bosonic environment.  
The dissipator is defined as follows
```math
D[J](\cdot) = J(\cdot) J^\dagger - \frac{1}{2}\left(J^\dagger J (\cdot) + (\cdot) J^\dagger J \right),
```
where ``J\equiv \sqrt{\gamma}V`` is the jump operator, ``V`` describes the dissipative part (operator) of the dynamics, ``\gamma`` represents a non-negative damping rate and ``[\cdot, \cdot]_+`` stands for anti-commutator.

Note that if ``V`` is acting on fermionic systems, it should be even-parity to be compatible with charge conservation.

# Parameters
- `M::AbstractHEOMLSMatrix` : the matrix given from HEOM model
- `jumpOP::AbstractVector` : The list of collapse (jump) operators ``\{J_i\}_i`` to add. Defaults to empty vector `[]`.

# Return 
- `M_new::AbstractHEOMLSMatrix` : the new HEOM Liouvillian superoperator matrix
"""
function addBosonDissipator(M::AbstractHEOMLSMatrix, jumpOP::Vector{T} = QuantumObject[]) where {T<:QuantumObject}
    if length(jumpOP) > 0
        return M + HEOMSuperOp(_sum_lindblad_dissipators(jumpOP), M.parity, M)
    else
        return M
    end
end
addBosonDissipator(M::AbstractHEOMLSMatrix, jumpOP::QuantumObject) = addBosonDissipator(M, [jumpOP])

@doc raw"""
    addFermionDissipator(M, jumpOP)

Adding fermionic dissipator to a given HEOMLS matrix which describes how the system dissipatively interacts with an extra fermionic environment.  
The dissipator with `EVEN` parity is defined as follows
```math
D_{\textrm{even}}[J](\cdot) = J(\cdot) J^\dagger - \frac{1}{2}\left(J^\dagger J (\cdot) + (\cdot) J^\dagger J \right),
```
where ``J\equiv \sqrt{\gamma}V`` is the jump operator, ``V`` describes the dissipative part (operator) of the dynamics, ``\gamma`` represents a non-negative damping rate and ``[\cdot, \cdot]_+`` stands for anti-commutator.

Similarly, the dissipator with `ODD` parity is defined as follows
```math
D_{\textrm{odd}}[J](\cdot) = - J(\cdot) J^\dagger - \frac{1}{2}\left(J^\dagger J (\cdot) + (\cdot) J^\dagger J \right),
```

Note that the parity of the dissipator will be determined by the parity of the given HEOMLS matrix `M`.

# Parameters
- `M::AbstractHEOMLSMatrix` : the matrix given from HEOM model
- `jumpOP::AbstractVector` : The list of collapse (jump) operators to add. Defaults to empty vector `[]`.

# Return 
- `M_new::AbstractHEOMLSMatrix` : the new HEOM Liouvillian superoperator matrix
"""
function addFermionDissipator(M::AbstractHEOMLSMatrix, jumpOP::Vector{T} = QuantumObject[]) where {T<:QuantumObject}
    if length(jumpOP) > 0
        L_data = mapreduce(J -> _fermion_lindblad_dissipator(J, M.parity), +, jumpOP)
        L = QuantumObject(L_data, type = SuperOperator(), dims = M.dimensions)

        return M + HEOMSuperOp(L, M.parity, M)
    else
        return M
    end
end
addFermionDissipator(M::AbstractHEOMLSMatrix, jumpOP::QuantumObject) = addFermionDissipator(M, [jumpOP])

function _fermion_lindblad_dissipator(J::QuantumObject{Operator}, parity::AbstractParity)
    _J = J.data
    Jd_J = _J' * _J
    return (-1)^(value(parity)) * _sprepost(_J, _J') - (_spre(Jd_J) + _spost(Jd_J)) / 2
end

@doc raw"""
    addTerminator(M, Bath)
Adding terminator to a given HEOMLS matrix.

The terminator is a Liouvillian term representing the contribution to 
the system-bath dynamics of all exponential-expansion terms beyond `Bath.Nterm`

The difference between the true correlation function and the sum of the 
`Bath.Nterm`-exponential terms is approximately `2 * δ * dirac(t)`.
Here, `δ` is the approximation discrepancy and `dirac(t)` denotes the Dirac-delta function.

# Parameters
- `M::AbstractHEOMLSMatrix` : the matrix given from HEOM model
- `Bath::Union{BosonBath, FermionBath}` : The bath object which contains the approximation discrepancy δ

# Return 
- `M_new::AbstractHEOMLSMatrix` : the new HEOM Liouvillian superoperator matrix
"""
function addTerminator(M::Mtype, Bath::Union{BosonBath,FermionBath}) where {Mtype<:AbstractHEOMLSMatrix}
    Btype = typeof(Bath)
    if (Btype == BosonBath) && (Mtype <: M_Fermion)
        error("For $(Btype), the type of HEOMLS matrix should be either M_Boson or M_Boson_Fermion.")
    elseif (Btype == FermionBath) && (Mtype <: M_Boson)
        error("For $(Btype), the type of HEOMLS matrix should be either M_Fermion or M_Boson_Fermion.")
    elseif Mtype <: M_S
        error("The type of input HEOMLS matrix does not support this functionality.")
    end

    if M.dimensions != Bath.op.dimensions
        error("The system dimensions between the HEOMLS matrix and Bath coupling operator are not consistent.")
    end

    if Bath.δ == 0
        @warn "The value of approximation discrepancy δ is 0.0 now, which doesn't make any changes."
        return M

    else
        return M + HEOMSuperOp(2 * Bath.δ * lindblad_dissipator(Bath.op), M.parity, M)
    end
end

raw"""
A sparse COO representation of the hierarchy (ADO) connectivity patterns and corresponding prefix values for a specific superoperator. The `data` field stores the corresponding superoperator matrix.
"""
struct COOFormat
    I::Vector{Int}
    J::Vector{Int}
    V::Vector{ComplexF64}
    N::Int
    data::AbstractMatrix
end

COOFormat(N::Int, data::AbstractMatrix) = COOFormat(Int[], Int[], ComplexF64[], N, data)

function Base.push!(mat::COOFormat, i, j, v)
    push!(mat.I, i)
    push!(mat.J, j)
    push!(mat.V, v)

    return nothing
end

SparseArrays.sparse(coo::COOFormat) = sparse(coo.I, coo.J, coo.V, coo.N, coo.N)

_gen_HEOMLS_term(coo::COOFormat) = TensorProductOperator(sparse(coo), coo.data)

raw"""
Stores the sparsity structure (positions and prefix values) of all superoperators in HEOM Liouville space using COO format.
"""
Base.@kwdef struct HEOMSparseStructure{
    Tspre<:Union{COOFormat,Nothing},
    Tspost<:Union{COOFormat,Nothing},
    TspreD<:Union{COOFormat,Nothing},
    TspostD<:Union{COOFormat,Nothing},
    TComm<:Union{COOFormat,Nothing},
    TanComm<:Union{COOFormat,Nothing},
    TCommD<:Union{COOFormat,Nothing},
}
    spre::Tspre = nothing
    spost::Tspost = nothing
    spreD::TspreD = nothing
    spostD::TspostD = nothing
    Comm::TComm = nothing
    anComm::TanComm = nothing
    CommD::TCommD = nothing
end

const HEOMSparseStructureFieldNames = fieldnames(HEOMSparseStructure)

HEOMSparseStructure(bath::AbstractFermionBath, Nado::Int) = HEOMSparseStructure(
    spre = COOFormat(Nado, bath.spre),
    spost = COOFormat(Nado, bath.spost),
    spreD = COOFormat(Nado, bath.spreD),
    spostD = COOFormat(Nado, bath.spostD),
)

HEOMSparseStructure(bath::bosonAbsorb, Nado::Int) = HEOMSparseStructure(
    spre = COOFormat(Nado, bath.spre),
    spost = COOFormat(Nado, bath.spost),
    CommD = COOFormat(Nado, bath.CommD),
)
HEOMSparseStructure(bath::bosonEmit, Nado::Int) = HEOMSparseStructure(
    spre = COOFormat(Nado, bath.spre),
    spost = COOFormat(Nado, bath.spost),
    CommD = COOFormat(Nado, bath.CommD),
)
HEOMSparseStructure(bath::bosonImag, Nado::Int) =
    HEOMSparseStructure(Comm = COOFormat(Nado, bath.Comm), anComm = COOFormat(Nado, bath.anComm))
HEOMSparseStructure(bath::bosonReal, Nado::Int) = HEOMSparseStructure(Comm = COOFormat(Nado, bath.Comm))
HEOMSparseStructure(bath::bosonRealImag, Nado::Int) =
    HEOMSparseStructure(Comm = COOFormat(Nado, bath.Comm), anComm = COOFormat(Nado, bath.anComm))

# sum γ of bath for current level
function bath_sum_γ(nvec, baths::Vector{T}) where {T<:Union{AbstractBosonBath,AbstractFermionBath}}
    p = 0
    sum_γ = 0.0
    for b in baths
        n = nvec[(p+1):(p+b.Nterm)]
        for k in findall(nk -> nk > 0, n)
            sum_γ += n[k] * b.γ[k]
        end
        p += b.Nterm
    end
    return sum_γ
end

# commutator of system Hamiltonian
minus_i_L_op(Hsys::QuantumObject) = liouvillian(Hsys).data

# connect to bosonic (n-1)th-level for "Real & Imag combined operator"
function minus_i_D_op!(ops_pattern::HEOMSparseStructure, I::Int, J::Int, bath::bosonRealImag, k, n_k)
    push!(ops_pattern.Comm, I, J, -1.0im * n_k * bath.η_real[k])
    push!(ops_pattern.anComm, I, J, n_k * bath.η_imag[k])
    return nothing
end

# connect to bosonic (n-1)th-level for (Real & Imag combined) operator "Real operator"
function minus_i_D_op!(ops_pattern::HEOMSparseStructure, I::Int, J::Int, bath::bosonReal, k, n_k)
    push!(ops_pattern.Comm, I, J, -1.0im * n_k * bath.η[k])
    return nothing
end

# connect to bosonic (n-1)th-level for "Imag operator"
function minus_i_D_op!(ops_pattern::HEOMSparseStructure, I::Int, J::Int, bath::bosonImag, k, n_k)
    push!(ops_pattern.anComm, I, J, n_k * bath.η[k])
    return nothing
end

# connect to bosonic (n-1)th-level for "Absorption operator"
function minus_i_D_op!(ops_pattern::HEOMSparseStructure, I::Int, J::Int, bath::bosonAbsorb, k, n_k)
    push!(ops_pattern.spre, I, J, -1.0im * n_k * bath.η[k])
    push!(ops_pattern.spost, I, J, 1.0im * n_k * conj(bath.η_emit[k]))
    return nothing
end

# connect to bosonic (n-1)th-level for "Emission operator"
function minus_i_D_op!(ops_pattern::HEOMSparseStructure, I::Int, J::Int, bath::bosonEmit, k, n_k)
    push!(ops_pattern.spre, I, J, -1.0im * n_k * bath.η[k])
    push!(ops_pattern.spost, I, J, 1.0im * n_k * conj(bath.η_absorb[k]))
    return nothing
end

# connect to fermionic (n-1)th-level for "absorption operator"
function minus_i_C_op!(
    ops_pattern::HEOMSparseStructure,
    I::Int,
    J::Int,
    bath::fermionAbsorb,
    k,
    n_exc,
    n_exc_before,
    parity,
)
    prefix = -1.0im * ((-1)^n_exc_before)
    push!(ops_pattern.spre, I, J, prefix * ((-1)^value(parity)) * bath.η[k])
    push!(ops_pattern.spost, I, J, - prefix * (-1)^(n_exc - 1) * conj(bath.η_emit[k]))
    return nothing
end

# connect to fermionic (n-1)th-level for "emission operator"
function minus_i_C_op!(
    ops_pattern::HEOMSparseStructure,
    I::Int,
    J::Int,
    bath::fermionEmit,
    k,
    n_exc,
    n_exc_before,
    parity,
)
    prefix = -1.0im * ((-1)^n_exc_before)
    push!(ops_pattern.spre, I, J, prefix * ((-1)^value(parity)) * bath.η[k])
    push!(ops_pattern.spost, I, J, - prefix * (-1)^(n_exc - 1) * conj(bath.η_absorb[k]))
    return nothing
end

# connect to bosonic (n+1)th-level for real-and-imaginary-type bosonic bath
function minus_i_B_op!(
    ops_pattern::HEOMSparseStructure,
    I::Int,
    J::Int,
    bath::T,
) where {T<:Union{bosonReal,bosonImag,bosonRealImag}}
    push!(ops_pattern.Comm, I, J, -1.0im)
    return nothing
end

# connect to bosonic (n+1)th-level for absorption-and-emission-type bosonic bath
function minus_i_B_op!(
    ops_pattern::HEOMSparseStructure,
    I::Int,
    J::Int,
    bath::T,
) where {T<:Union{bosonAbsorb,bosonEmit}}
    push!(ops_pattern.CommD, I, J, -1.0im)
    return nothing
end

# connect to fermionic (n+1)th-level
function minus_i_A_op!(
    ops_pattern::HEOMSparseStructure,
    I::Int,
    J::Int,
    bath::AbstractFermionBath,
    n_exc,
    n_exc_before,
    parity,
)
    prefix = -1.0im * ((-1)^n_exc_before)
    push!(ops_pattern.spreD, I, J, prefix * ((-1)^value(parity)))
    push!(ops_pattern.spostD, I, J, prefix * (-1)^(n_exc + 1))
    return nothing
end

function combine_HEOMLS_terms(op::AddedOperator)
    Tensor_ops = op.ops |> collect # [ A_i ⊗ B_i ]
    A_list = Vector{AbstractMatrix}(undef, length(Tensor_ops))
    B_list = Vector{AbstractMatrix}(undef, length(Tensor_ops))
    for (idx, t_op) in pairs(Tensor_ops)
        t_op isa TensorProductOperator || throw(ArgumentError("The HEOMLS term should be a TensorProductOperator."))
        A_list[idx] = t_op.ops[1].A
        B_list[idx] = t_op.ops[2].A
    end

    unique_B_ops = unique(B_list)
    index_groups = [[] for i in unique_B_ops]
    for i in eachindex(Tensor_ops)
        for j in eachindex(unique_B_ops)
            if isequal(unique_B_ops[j], B_list[i])
                push!(index_groups[j], i)
            end
        end
    end

    return sum(pairs(unique_B_ops)) do (j, Bj)
        Aj = sum(k -> A_list[k], index_groups[j])
        return TensorProductOperator(Aj, Bj)
    end
end
combine_HEOMLS_terms(op::TensorProductOperator) = op

function assemble_HEOMLS_terms(M::Vector{<:AbstractSciMLOperator}, ::Val{:full}, verbose::Bool)
    M_combine = assemble_HEOMLS_terms(M, Val(:combine), verbose) # combine first

    if verbose
        print("Evaluating lazy operations...")
        flush(stdout)
    end
    M_full = map(x -> MatrixOperator(SciMLOperators.concretize(x)), M_combine)
    if verbose
        println("[DONE]")
        flush(stdout)
    end
    return M_full
end
function assemble_HEOMLS_terms(M::Vector{<:AbstractSciMLOperator}, ::Val{:combine}, verbose::Bool)
    if verbose
        print("Combining terms...")
        flush(stdout)
    end
    M_combine = map(combine_HEOMLS_terms, M)
    if verbose
        println("[DONE]")
        flush(stdout)
    end
    return M_combine
end
assemble_HEOMLS_terms(M::Vector{<:AbstractSciMLOperator}, ::Val{:none}, verbose::Bool) = M
assemble_HEOMLS_terms(M::AbstractSciMLOperator, method::Val, verbose::Bool) =
    assemble_HEOMLS_terms([M], method, verbose)

check_assemble_method(assemble_method) =
    (assemble_method ∉ (Val(:full), Val(:combine), Val(:none))) && throw(
        ArgumentError(
            "Invalid value for `assemble`: $(assemble_method). Accepted values are `:full`, `:combine`, or `:none`.",
        ),
    )

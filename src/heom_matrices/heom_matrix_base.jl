export HEOMSuperOp
export Propagator
export addBosonDissipator, addFermionDissipator, addTerminator

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
    HEOMSuperOp(op, opParity, refHEOMLS; Id_cache=I(refHEOMLS.N))
Construct the HEOM superoperator matrix corresponding to the given system SuperOperator which acts on all `ADOs`.  

During the multiplication on all the `ADOs`, the parity of the output `ADOs` might change depend on the parity of this HEOM superoperator.

# Parameters
- `op` : The system SuperOperator which will act on all `ADOs`.
- `opParity::AbstractParity` : the parity label of the given operator (`op`), should be `EVEN` or `ODD`.
- `refHEOMLS::AbstractHEOMLSMatrix` : copy the system `dimensions` and number of `ADOs` (`N`) from this reference HEOMLS matrix
"""
HEOMSuperOp(op, opParity::AbstractParity, refHEOMLS::AbstractHEOMLSMatrix; Id_cache = I(refHEOMLS.N)) =
    HEOMSuperOp(op, opParity, refHEOMLS.dimensions, refHEOMLS.N; Id_cache = Id_cache)

@doc raw"""
    HEOMSuperOp(op, opParity, refADOs; Id_cache=I(refADOs.N))
Construct the HEOM superoperator matrix corresponding to the given system SuperOperator which acts on all `ADOs`.  

During the multiplication on all the `ADOs`, the parity of the output `ADOs` might change depend on the parity of this HEOM superoperator.

# Parameters
- `op` : The system SuperOperator which will act on all `ADOs`.
- `opParity::AbstractParity` : the parity label of the given operator (`op`), should be `EVEN` or `ODD`.
- `refADOs::ADOs` : copy the system `dimensions` and number of `ADOs` (`N`) from this reference `ADOs`   
"""
HEOMSuperOp(op, opParity::AbstractParity, refADOs::ADOs; Id_cache = I(refADOs.N)) =
    HEOMSuperOp(op, opParity, refADOs.dimensions, refADOs.N; Id_cache = Id_cache)

@doc raw"""
    HEOMSuperOp(op, opParity, dims, N; Id_cache=I(N))
Construct the HEOM superoperator matrix corresponding to the given system SuperOperator which acts on all `ADOs`.  

During the multiplication on all the `ADOs`, the parity of the output `ADOs` might change depend on the parity of this HEOM superoperator.

# Parameters
- `op` : The system SuperOperator which will act on all `ADOs`.
- `opParity::AbstractParity` : the parity label of the given operator (`op`), should be `EVEN` or `ODD`.
- `dims` : the dimension list of the coupling operator (should be equal to the system `dimensions`).
- `N::Int` : the number of `ADOs`.
"""
function HEOMSuperOp(op, opParity::AbstractParity, dims, N::Int; Id_cache = I(N))
    dimensions = _gen_dimensions(dims)
    sup_op = HandleMatrixType(op, dimensions, "op (operator)"; type = SuperOperator())

    return HEOMSuperOp(kron(Id_cache, sup_op.data), dimensions, N, opParity)
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
        Id_cache = I(prod(M.dimensions))
        L = mapreduce(J -> lindblad_dissipator(J, Id_cache), +, jumpOP)

        return M + HEOMSuperOp(L, M.parity, M)
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
        Id_cache = I(prod(M.dimensions))
        L_data = mapreduce(J -> _fermion_lindblad_dissipator(J, M.parity, Id_cache), +, jumpOP)
        L = QuantumObject(L_data, type = SuperOperator(), dims = M.dimensions)

        return M + HEOMSuperOp(L, M.parity, M)
    else
        return M
    end
end
addFermionDissipator(M::AbstractHEOMLSMatrix, jumpOP::QuantumObject) = addFermionDissipator(M, [jumpOP])

function _fermion_lindblad_dissipator(J::QuantumObject{Operator}, parity::AbstractParity, Id_cache = I(size(J, 1)))
    _J = J.data
    Jd_J = _J' * _J
    return (-1)^(value(parity)) * _sprepost(_J, _J') - (_spre(Jd_J, Id_cache) + _spost(Jd_J, Id_cache)) / 2
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

function csc2coo(A)
    len = length(A.nzval)

    if len == 0
        return A.m, A.n, [], [], []
    else
        colidx = Vector{Int}(undef, len)
        @inbounds for i in 1:(length(A.colptr)-1)
            @inbounds for j in A.colptr[i]:(A.colptr[i+1]-1)
                colidx[j] = i
            end
        end
        return A.m, A.n, A.rowval, colidx, A.nzval
    end
end

function pad_coo(
    A::SparseMatrixCSC{T,Int64},
    row_scale::Int,
    col_scale::Int,
    row_idx = 1::Int,
    col_idx = 1::Int,
) where {T<:Number}
    # transform matrix A's format from csc to coo
    M, N, I, J, V = csc2coo(A)

    # deal with values
    if T != ComplexF64
        V = convert.(ComplexF64, V)
    end

    # deal with rowval
    if (row_idx > row_scale) || (row_idx < 1)
        error("row_idx must be \'>= 1\' and \'<= row_scale\'")
    end

    # deal with colval
    if (col_idx > col_scale) || (col_idx < 1)
        error("col_idx must be \'>= 1\' and \'<= col_scale\'")
    end

    @inbounds Inew = I .+ (M * (row_idx - 1))
    @inbounds Jnew = J .+ (N * (col_idx - 1))

    return Inew, Jnew, V
end

function add_operator!(op, I, J, V, N_he, row_idx, col_idx)
    row, col, val = pad_coo(op, N_he, N_he, row_idx, col_idx)
    append!(I, row)
    append!(J, col)
    append!(V, val)
    return nothing
end

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
minus_i_L_op(Hsys::QuantumObject, Id = I(size(Hsys, 1))) = liouvillian(Hsys, Id).data

# connect to bosonic (n-1)th-level for "Real & Imag combined operator"
minus_i_D_op(bath::bosonRealImag, k, n_k) = n_k * (-1.0im * bath.η_real[k] * bath.Comm + bath.η_imag[k] * bath.anComm)

# connect to bosonic (n-1)th-level for (Real & Imag combined) operator "Real operator"
minus_i_D_op(bath::bosonReal, k, n_k) = -1.0im * n_k * bath.η[k] * bath.Comm

# connect to bosonic (n-1)th-level for "Imag operator"
minus_i_D_op(bath::bosonImag, k, n_k) = n_k * bath.η[k] * bath.anComm

# connect to bosonic (n-1)th-level for "Absorption operator"
minus_i_D_op(bath::bosonAbsorb, k, n_k) = -1.0im * n_k * (bath.η[k] * bath.spre - conj(bath.η_emit[k]) * bath.spost)

# connect to bosonic (n-1)th-level for "Emission operator"
minus_i_D_op(bath::bosonEmit, k, n_k) = -1.0im * n_k * (bath.η[k] * bath.spre - conj(bath.η_absorb[k]) * bath.spost)

# connect to fermionic (n-1)th-level for "absorption operator"
function minus_i_C_op(bath::fermionAbsorb, k, n_exc, n_exc_before, parity)
    return -1.0im *
           ((-1)^n_exc_before) *
           (((-1)^value(parity)) * bath.η[k] * bath.spre - (-1)^(n_exc - 1) * conj(bath.η_emit[k]) * bath.spost)
end

# connect to fermionic (n-1)th-level for "emission operator"
function minus_i_C_op(bath::fermionEmit, k, n_exc, n_exc_before, parity)
    return -1.0im *
           ((-1)^n_exc_before) *
           ((-1)^(value(parity)) * bath.η[k] * bath.spre - (-1)^(n_exc - 1) * conj(bath.η_absorb[k]) * bath.spost)
end

# connect to bosonic (n+1)th-level for real-and-imaginary-type bosonic bath
minus_i_B_op(bath::T) where {T<:Union{bosonReal,bosonImag,bosonRealImag}} = -1.0im * bath.Comm

# connect to bosonic (n+1)th-level for absorption-and-emission-type bosonic bath
minus_i_B_op(bath::T) where {T<:Union{bosonAbsorb,bosonEmit}} = -1.0im * bath.CommD

# connect to fermionic (n+1)th-level
function minus_i_A_op(bath::T, n_exc, n_exc_before, parity) where {T<:AbstractFermionBath}
    return -1.0im * ((-1)^n_exc_before) * ((-1)^(value(parity)) * bath.spreD + (-1)^(n_exc + 1) * bath.spostD)
end

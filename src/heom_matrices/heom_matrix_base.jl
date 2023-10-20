abstract type AbstractHEOMLSMatrix end

# Parity label
const odd  = 1;
const even = 0;

@doc raw"""
    size(M::AbstractHEOMLSMatrix)
Returns the size of the HEOM Liouvillian superoperator matrix
"""
size(M::AbstractHEOMLSMatrix) = size(M.data)

@doc raw"""
    size(M::AbstractHEOMLSMatrix, dim::Int)
Returns the specified dimension of the HEOM Liouvillian superoperator matrix
"""
size(M::AbstractHEOMLSMatrix, dim::Int) = size(M.data, dim)

getindex(M::AbstractHEOMLSMatrix, i::Ti, j::Tj) where {Ti, Tj <: Any} = M.data[i, j]

function show(io::IO, M::AbstractHEOMLSMatrix)
    T = typeof(M)
    if T == M_S
        type = "Schrodinger Eq."
    elseif T == M_Boson
        type = "Boson"
    elseif T == M_Fermion
        type = "Fermion"
    else
        type = "Boson-Fermion"
    end

    print(io, 
        type, " type HEOMLS matrix with (system) dim = $(M.dim) and parity = :$(M.parity)\n",
        "number of ADOs N = $(M.N)\n",
        "data =\n"
    )
    show(io, MIME("text/plain"), M.data)
end

function show(io::IO, m::MIME"text/plain", M::AbstractHEOMLSMatrix) show(io, M) end

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
        M::AbstractHEOMLSMatrix,
        Δt::Real;
        threshold   = 1.0e-6,
        nonzero_tol = 1.0e-14
    )
    return fastExpm(M.data * Δt; threshold=threshold, nonzero_tol=nonzero_tol)
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
function addBosonDissipator(M::T, jumpOP::Vector=[]) where T <: AbstractHEOMLSMatrix
    if length(jumpOP) > 0
        L = spzeros(ComplexF64, M.sup_dim, M.sup_dim)
        for J in jumpOP
            _J = HandleMatrixType(J, M.dim, "in jumpOP")
            L += spre(_J) * spost(_J') - 0.5 * (spre(_J' * _J) + spost(_J' * _J))
        end

        if T == M_S
            return M_S(M.data + L, M.tier, M.dim, M.N, M.sup_dim, M.parity)
        elseif T == M_Boson
            return M_Boson(M.data + kron(sparse(I, M.N, M.N), L), M.tier, M.dim, M.N, M.sup_dim, M.parity, M.bath, M.hierarchy)
        elseif T == M_Fermion
            return M_Fermion(M.data + kron(sparse(I, M.N, M.N), L), M.tier, M.dim, M.N, M.sup_dim, M.parity, M.bath, M.hierarchy)
        else
            return M_Boson_Fermion(M.data + kron(sparse(I, M.N, M.N), L), M.Btier, M.Ftier, M.dim, M.N, M.sup_dim, M.parity, M.Bbath, M.Fbath, M.hierarchy)
        end
    end
end

@doc raw"""
    addFermionDissipator(M, jumpOP)
Adding fermionic dissipator to a given HEOMLS matrix which describes how the system dissipatively interacts with an extra fermionic environment.  
The dissipator with `:even` parity is defined as follows
```math
D_{\textrm{even}}[J](\cdot) = J(\cdot) J^\dagger - \frac{1}{2}\left(J^\dagger J (\cdot) + (\cdot) J^\dagger J \right),
```
where ``J\equiv \sqrt{\gamma}V`` is the jump operator, ``V`` describes the dissipative part (operator) of the dynamics, ``\gamma`` represents a non-negative damping rate and ``[\cdot, \cdot]_+`` stands for anti-commutator.

Similary, the dissipator with `:odd` parity is defined as follows
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
function addFermionDissipator(M::T, jumpOP::Vector=[]) where T <: AbstractHEOMLSMatrix
    if length(jumpOP) > 0
        parity = eval(M.parity)
        L = spzeros(ComplexF64, M.sup_dim, M.sup_dim)
        for J in jumpOP
            _J = HandleMatrixType(J, M.dim, "in jumpOP")
            L += ((-1) ^ parity) * spre(_J) * spost(_J') - 0.5 * (spre(_J' * _J) + spost(_J' * _J))
        end
        if T == M_S
            return M_S(M.data + L, M.tier, M.dim, M.N, M.sup_dim, M.parity)
        elseif T == M_Boson
            return M_Boson(M.data + kron(sparse(I, M.N, M.N), L), M.tier, M.dim, M.N, M.sup_dim, M.parity, M.bath, M.hierarchy)
        elseif T == M_Fermion
            return M_Fermion(M.data + kron(sparse(I, M.N, M.N), L), M.tier, M.dim, M.N, M.sup_dim, M.parity, M.bath, M.hierarchy)
        else
            return M_Boson_Fermion(M.data + kron(sparse(I, M.N, M.N), L), M.Btier, M.Ftier, M.dim, M.N, M.sup_dim, M.parity, M.Bbath, M.Fbath, M.hierarchy)
        end
    end
end

function   addBosonDissipator(M::AbstractHEOMLSMatrix, jumpOP) return   addBosonDissipator(M, [jumpOP]) end
function addFermionDissipator(M::AbstractHEOMLSMatrix, jumpOP) return addFermionDissipator(M, [jumpOP]) end

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
function addTerminator(M::Mtype, Bath::Union{BosonBath, FermionBath}) where Mtype <: AbstractHEOMLSMatrix
    Btype = typeof(Bath)
    if (Btype == BosonBath) && (Mtype == M_Fermion)
        error("For $(Btype), the type of HEOMLS matrix should be either M_Boson or M_Boson_Fermion.")
    elseif (Btype == FermionBath) && (Mtype == M_Boson)
        error("For $(Btype), the type of HEOMLS matrix should be either M_Fermion or M_Boson_Fermion.")
    elseif Mtype == M_S
        error("The type of input HEOMLS matrix does not support this functionality.")
    end

    if M.dim != Bath.dim
        error("The system dimension between the HEOMLS matrix and Bath are not consistent.")
    end

    if Bath.δ == 0
        @warn "The value of approximation discrepancy δ is 0.0 now, which doesn't make any changes."
        return M

    else
        J = Bath.op
        L = 2 * Bath.δ * (
            spre(J) * spost(J') - 0.5 * (spre(J' * J) + spost(J' * J))
        )

        if Mtype == M_Boson
            return M_Boson(M.data + kron(sparse(I, M.N, M.N), L), M.tier, M.dim, M.N, M.sup_dim, M.parity, M.bath, M.hierarchy)
        elseif Mtype == M_Fermion
            return M_Fermion(M.data + kron(sparse(I, M.N, M.N), L), M.tier, M.dim, M.N, M.sup_dim, M.parity, M.bath, M.hierarchy)
        else
            return M_Boson_Fermion(M.data + kron(sparse(I, M.N, M.N), L), M.Btier, M.Ftier, M.dim, M.N, M.sup_dim, M.parity, M.Bbath, M.Fbath, M.hierarchy)
        end
    end
end

function csc2coo(A)
    len = length(A.nzval)

    if len == 0
        return A.m, A.n, [], [], []
    else        
        colidx = Vector{Int}(undef, len)
        @inbounds for i in 1:(length(A.colptr) - 1)
            @inbounds for j in A.colptr[i] : (A.colptr[i + 1] - 1)
                colidx[j] = i
            end
        end
        return A.m, A.n, A.rowval, colidx, A.nzval
    end
end

function pad_coo(A::SparseMatrixCSC{T, Int64}, row_scale::Int, col_scale::Int, row_idx=1::Int, col_idx=1::Int) where {T<:Number}
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
    nothing
end

# sum γ of bath for current level
function bath_sum_γ(nvec, baths::Vector{T}) where T <: Union{AbstractBosonBath, AbstractFermionBath}
    p = 0
    sum_γ = 0.0
    for b in baths
        n = nvec[(p + 1) : (p + b.Nterm)]
        for k in findall(nk -> nk > 0, n)
            sum_γ += n[k] * b.γ[k]
        end
        p += b.Nterm
    end
    return sum_γ
end

# connect to bosonic (n-1)th-level for "Real & Imag combined operator"
function _D_op(bath::bosonRealImag, k, n_k)
    return n_k * (
        -1im * bath.η_real[k] * bath.Comm  +
               bath.η_imag[k] * bath.anComm
    )
end

# connect to bosonic (n-1)th-level for (Real & Imag combined) operator "Real operator"
function _D_op(bath::bosonReal, k, n_k)
    return -1im * n_k * bath.η[k] * bath.Comm
end

# connect to bosonic (n-1)th-level for "Imag operator"
function _D_op(bath::bosonImag, k, n_k)
    return n_k * bath.η[k] * bath.anComm
end

# connect to bosonic (n-1)th-level for "Absorption operator"
function _D_op(bath::bosonAbsorb, k, n_k)
    return -1im * n_k * (
        bath.η[k] * bath.spre - conj(bath.η_emit[k]) * bath.spost
    )
end

# connect to bosonic (n-1)th-level for "Emission operator"
function _D_op(bath::bosonEmit, k, n_k)
    return -1im * n_k * (
        bath.η[k] * bath.spre - conj(bath.η_absorb[k]) * bath.spost
    )
end

# connect to fermionic (n-1)th-level for "absorption operator"
function _C_op(bath::fermionAbsorb, k, n_exc, n_exc_before, parity)
    return -1im * ((-1) ^ n_exc_before) * (
        ((-1) ^ eval(parity)) * bath.η[k] * bath.spre -
        (-1) ^ (n_exc - 1)    * conj(bath.η_emit[k]) * bath.spost
    )
end

# connect to fermionic (n-1)th-level for "emission operator"
function _C_op(bath::fermionEmit, k, n_exc, n_exc_before, parity)
    return -1im * ((-1) ^ n_exc_before) * (
        (-1) ^ (eval(parity)) * bath.η[k] * bath.spre -
        (-1) ^ (n_exc - 1)    * conj(bath.η_absorb[k]) * bath.spost
    )
end

# connect to bosonic (n+1)th-level for real-and-imaginary-type bosonic bath
function _B_op(bath::T) where T <: Union{bosonReal, bosonImag, bosonRealImag}
    return -1im * bath.Comm
end

# connect to bosonic (n+1)th-level for absorption-and-emission-type bosonic bath
function _B_op(bath::T) where T <: Union{bosonAbsorb, bosonEmit}
    return -1im * bath.CommD
end

# connect to fermionic (n+1)th-level
function _A_op(bath::T, n_exc, n_exc_before, parity) where T <: AbstractFermionBath
    return -1im * ((-1) ^ n_exc_before) * (
        (-1) ^ (eval(parity)) * bath.spreD  +
        (-1) ^ (n_exc + 1)    * bath.spostD
    )
end
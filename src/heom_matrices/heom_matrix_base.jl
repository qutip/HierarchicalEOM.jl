abstract type AbstractHEOMMatrix end

# Parity for fermionic heom matrices
const odd  = 1;
const even = 0;
const none = nothing;

"""
    size(M::AbstractHEOMMatrix)
Returns the size of the Heom liouvillian superoperator matrix
"""
size(M::AbstractHEOMMatrix) = size(M.data)

getindex(M::AbstractHEOMMatrix, i::Ti, j::Tj) where {Ti, Tj <: Any} = M.data[i, j]

function show(io::IO, M::AbstractHEOMMatrix)
    T = typeof(M)
    if T == M_Boson
        type = "Boson"
    elseif T == M_Fermion
        type = "Fermion"
    else
        type = "Boson-Fermion"
    end

    print(io, 
        type, " type HEOM matrix with (system) dim = $(M.dim) and parity = :$(M.parity)\n",
        "number of ADOs N = $(M.N)\n",
        "data =\n"
    )
    show(io, MIME("text/plain"), M.data)
end

function show(io::IO, m::MIME"text/plain", M::AbstractHEOMMatrix) show(io, M) end

"""
    Propagator(M, Δt; threshold, nonzero_tol)
Use `FastExpm.jl` to calculate the propagator matrix from a given Heom liouvillian superoperator matrix ``M`` with a specific time step ``Δt``.
That is, ``\\exp(M * \\Delta t)``.

# Parameters
- `M::AbstractHEOMMatrix` : the matrix given from HEOM model
- `Δt::Real` : A specific time step (time interval).
- `threshold::Real` : Determines the threshold for the Taylor series. Defaults to `1.0e-6`.
- `nonzero_tol::Real` : Strips elements smaller than `nonzero_tol` at each computation step to preserve sparsity. Defaults to `1.0e-14`.

For more details, please refer to [`FastExpm.jl`](https://github.com/fmentink/FastExpm.jl)

# Returns
- `::SparseMatrixCSC{ComplexF64, Int64}` : the propagator matrix
"""
function Propagator(
        M::AbstractHEOMMatrix,
        Δt::Real;
        threshold   = 1.0e-6,
        nonzero_tol = 1.0e-14
    )
    return fastExpm(M.data * Δt; threshold=threshold, nonzero_tol=nonzero_tol)
end

"""
    addDissipator!(M, jumpOP)
Adding dissipator to a given HEOM matrix.

# Parameters
- `M::AbstractHEOMMatrix` : the matrix given from HEOM model
- `jumpOP::AbstractVector` : The collapse (jump) operators to add. Defaults to empty vector `[]`.
"""
function addDissipator!(M::AbstractHEOMMatrix, jumpOP::Vector=[])
    if length(jumpOP) > 0
        L = spzeros(ComplexF64, M.sup_dim, M.sup_dim)
        for J in jumpOP
            if isValidMatrixType(J, M.dim)
                L += spre(J) * spost(J') - 0.5 * (spre(J' * J) + spost(J' * J))
            else
                error("Invalid matrix in \"jumpOP\".")
            end
        end
        M.data += kron(sparse(I, M.N, M.N), L)
    end
end

function addDissipator!(M::AbstractHEOMMatrix, jumpOP) addDissipator!(M, [jumpOP]) end

"""
    addTerminator!(M, Bath)
Adding terminator to a given HEOM matrix.

The terminator is a liouvillian term representing the contribution to 
the system-bath dynamics of all exponential-expansion terms beyond `Bath.Nterm`

The difference between the true correlation function and the sum of the 
`Bath.Nterm`-exponential terms is approximately `2 * δ * dirac(t)`.
Here, `δ` is the approximation discrepancy and `dirac(t)` denotes the Dirac-delta function.

# Parameters
- `M::AbstractHEOMMatrix` : the matrix given from HEOM model
- `Bath::Union{BosonBath, FermionBath}` : The bath object which contains the approximation discrepancy δ
"""
function addTerminator!(M::AbstractHEOMMatrix, Bath::Union{BosonBath, FermionBath})
    Btype = typeof(Bath)
    Mtype = typeof(M)
    if (Btype == BosonBath) && (Mtype == M_Fermion)
        error("For $(Btype), the type of Heom matrix should be either M_Boson or M_Boson_Fermion.")
    elseif (Btype == FermionBath) && (Mtype == M_Boson)
        error("For $(Btype), the type of Heom matrix should be either M_Fermion or M_Boson_Fermion.")
    end

    if M.dim != Bath.dim
        error("The system dimension between the Heom matrix and Bath are not consistent.")
    end

    if Bath.δ == 0
        @warn "The value of approximation discrepancy δ is 0.0 now, which doesn't make any changes."
    else
        J = Bath.op
        L = 2 * Bath.δ * (
            spre(J) * spost(J') - 0.5 * (spre(J' * J) + spost(J' * J))
        )

        M.data += kron(sparse(I, M.N, M.N), L)
    end
end

function pad_csc(A::SparseMatrixCSC{T, Int64}, row_scale::Int, col_scale::Int, row_idx=1::Int, col_idx=1::Int) where {T<:Number}
    (M, N) = size(A)

    # deal with values
    values = A.nzval
    if length(values) == 0
        return sparse([M * row_scale], [N * col_scale], [0.0im])
    else
        if T != ComplexF64
            values = convert.(ComplexF64, values)
        end

        # deal with colptr
        local ptrLen::Int         = N * col_scale + 1
        local ptrIn::Vector{Int}  = A.colptr
        local ptrOut::Vector{Int} = fill(1, ptrLen)
        if col_idx == 1
            ptrOut[1:(N+1)]   .= ptrIn            
            ptrOut[(N+2):end] .= ptrIn[end]

        elseif col_idx == col_scale         
            ptrOut[(ptrLen-N):end] .= ptrIn

        elseif (col_idx < col_scale) && (col_idx > 1)
            tmp1 = (col_idx - 1) * N + 1
            tmp2 = tmp1 + N
            ptrOut[tmp1:tmp2] .= ptrIn
            ptrOut[(tmp2+1):end] .= ptrIn[end]

        else
            error("col_idx must be \'>= 1\' and \'<= col_scale\'")
        end

        # deal with rowval
        if (row_idx > row_scale) || (row_idx < 1)
            error("row_idx must be \'>= 1\' and \'<= row_scale\'")
        end
        tmp1 = (row_idx - 1) * N

        return SparseMatrixCSC(
            M * row_scale,
            N * col_scale,
            ptrOut,
            A.rowval .+ tmp1, 
            values,
        )
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
    append!(localpart(I)[1], row)
    append!(localpart(J)[1], col)
    append!(localpart(V)[1], val)
end

# sum ω of bath for current gradient
function bath_sum_ω(nvec, baths::Vector{T}) where T <: Union{AbstractBosonBath, AbstractFermionBath}
    p = 0
    sum_ω = 0.0
    for b in baths
        for k in findall(nk -> nk > 0, nvec[(p + 1) : (p + b.Nterm)])
            sum_ω += nvec[p + k] * b.γ[k]
        end
        p += b.Nterm
    end
    return sum_ω
end

# boson (Real & Imag combined) operator for previous gradient
function prev_grad_boson(bath::bosonRealImag, k, n_k)
    return n_k * (
        -1im * bath.η_real[k] * bath.Comm  +
               bath.η_imag[k] * bath.anComm
    )
end

# boson (Real) operator for previous gradient
function prev_grad_boson(bath::bosonReal, k, n_k)
    return -1im * n_k * bath.η[k] * bath.Comm
end

# boson (Imag) operator for previous gradient
function prev_grad_boson(bath::bosonImag, k, n_k)
    return n_k * bath.η[k] * bath.anComm
end

# absorption fermion operator for previous gradient
function prev_grad_fermion(bath::fermionAbsorb, k, n_exc, n_exc_before, parity)
    return -1im * ((-1) ^ n_exc_before) * (
        (-1) ^ eval(parity) * bath.η[k] * bath.spre - 
        (-1) ^ (n_exc - 1)  * conj(bath.η_emit[k]) * bath.spost
    )
end

# emission fermion operator for previous gradient
function prev_grad_fermion(bath::fermionEmit, k, n_exc, n_exc_before, parity)
    return -1im * ((-1) ^ n_exc_before) * (
        (-1) ^ eval(parity) * bath.η[k] * bath.spre - 
        (-1) ^ (n_exc - 1)  * conj(bath.η_absorb[k]) * bath.spost
    )
end

# boson operator for next gradient
function next_grad_boson(bath::T) where T <: AbstractBosonBath
    return -1im * bath.Comm
end

# fermion operator for next gradient
function next_grad_fermion(bath::T, n_exc, n_exc_before, parity) where T <: AbstractFermionBath
    return -1im * ((-1) ^ n_exc_before) * (
        (-1) ^ eval(parity) * bath.spreD +
        (-1) ^ (n_exc - 1)  * bath.spostD
    )
end
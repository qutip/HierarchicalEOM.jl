abstract type AbstractHEOMMatrix end

# Parity for fermionic heom matrices
const odd  = 1;
const even = 0;
const none = nothing;

size(A::AbstractHEOMMatrix) = size(A.data)

spre(q::AbstractMatrix)        = kron(Matrix(I, size(q)[1], size(q)[1]), q)
spre(q::AbstractSparseMatrix)  = sparse(kron(sparse(I, size(q)[1], size(q)[1]), q))
spost(q::AbstractMatrix)       = kron(transpose(q), Matrix(I, size(q)[1], size(q)[1]))
spost(q::AbstractSparseMatrix) = sparse(kron(transpose(q), sparse(I, size(q)[1], size(q)[1])))

"""
# `addDissipator!(M, jumpOP)`
Adding dissipator to a given HEOM matrix.

## Parameters
- `M::AbstractHEOMMatrix` : the matrix given from HEOM model
- `jumpOP::Vector{T<:AbstractMatrix}` : The collapse (jump) operators to add. Defaults to empty vector `[]`.
"""
function addDissipator!(M::AbstractHEOMMatrix, jumpOP::Vector{T}=[]) where T <: AbstractMatrix
    if length(jumpOP) > 0
        for J in jumpOP
            M.data += spre(J) * spost(J') - 0.5 * (spre(J' * J) + spost(J' * J))
        end
    end
end

# generate index to ado vector
function ADO_number(dims::Vector{Int}, N_exc::Int)
    len = length(dims)
    state = zeros(Int, len)
    result = [copy(state)]
    nexc = 0

    while true
        idx = len
        state[end] += 1
        nexc += 1
        if state[idx] < dims[idx]
            push!(result, copy(state))
        end
        while (nexc == N_exc) || (state[idx] == dims[idx])
            #state[idx] = 0
            idx -= 1
            if idx < 1
                return result
            end

            nexc -= state[idx + 1] - 1
            state[idx + 1] = 0
            state[idx] += 1
            if state[idx] < dims[idx]
                push!(result, copy(state))
            end
        end
    end
end

function ADOs_dictionary(dims::Vector{Int}, N_exc::Int)
    ado2idx = OrderedDict{Vector{Int}, Int}()
    idx2ado = ADO_number(dims, N_exc)
    for (idx, ado) in enumerate(idx2ado)
        ado2idx[ado] = idx
    end

    return length(idx2ado), ado2idx, idx2ado
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
        colptr = A.colptr
        colidx = Vector{Int}(undef, len)

        @inbounds for i in 1:(length(A.colptr) - 1)
            @inbounds for j in colptr[i] : (colptr[i + 1] - 1)
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

    @inbounds I .+= (M * (row_idx - 1))
    @inbounds J .+= (N * (col_idx - 1))
    
    return I, J, V
end
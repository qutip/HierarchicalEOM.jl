module HierarchicalEOM_CUDAExt

using HierarchicalEOM
import HierarchicalEOM: _HandleVectorType, _HandleTraceVectorType
import QuantumToolbox: _CType, _convert_eltype_wordsize, makeVal, getVal
import CUDA
import CUDA: cu, CuArray
import CUDA.CUSPARSE: CuSparseVector, CuSparseMatrixCSC
import SparseArrays: AbstractSparseMatrix, sparse, SparseVector, SparseMatrixCSC
import SciMLOperators: MatrixOperator, ScaledOperator, AddedOperator

@doc raw"""
    cu(M::AbstractHEOMLSMatrix; word_size::Union{Val,Int} = Val(64))
Return a new HEOMLS-matrix-type object with `M.data` is in the type of `CuSparseMatrixCSC{ComplexF32, Int32}` for gpu calculations.

# Arguments
- `M::AbstractHEOMLSMatrix` : the matrix given from HEOM model
- `word_size::Union{Val,Int}` : The word size of the element type of `M`, can be either `32` or `64`. Default to `64`.
"""
function cu(M::AbstractHEOMLSMatrix; word_size::Union{Val,Int} = Val(64))
    word_size_val = makeVal(word_size)
    _word_size = getVal(word_size_val)

    ((_word_size == 64) || (_word_size == 32)) || throw(DomainError(_word_size, "The word size should be 32 or 64."))

    return CuSparseMatrixCSC{_convert_eltype_wordsize(eltype(M), word_size_val)}(M)
end

@doc raw"""
    CuSparseMatrixCSC{T}(M::AbstractHEOMLSMatrix)

Return a new HEOMLS-matrix-type object with `M.data` is in the type of `CUDA.CUSPARSE.CuSparseMatrixCSC` with element type `T` for gpu calculations.
"""
function CuSparseMatrixCSC{T}(M::MT) where {T,MT<:AbstractHEOMLSMatrix}
    A_gpu = _convert_to_gpu_matrix(M.data, T)
    if MT <: M_S
        return M_S(A_gpu, M.tier, M.dimensions, M.N, M.sup_dim, M.parity)
    elseif MT <: M_Boson
        return M_Boson(A_gpu, M.tier, M.dimensions, M.N, M.sup_dim, M.parity, M.bath, M.hierarchy)
    elseif MT <: M_Fermion
        return M_Fermion(A_gpu, M.tier, M.dimensions, M.N, M.sup_dim, M.parity, M.bath, M.hierarchy)
    else
        return M_Boson_Fermion(
            A_gpu,
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

CuSparseMatrixCSC{T}(M::HEOMSuperOp) where {T} =
    HEOMSuperOp(_convert_to_gpu_matrix(M.data, T), M.dimensions, M.N, M.parity)

function _convert_to_gpu_matrix(A::AbstractSparseMatrix, ElType::Type{T}) where {T<:Number}
    if A isa CuSparseMatrixCSC{ElType,Int32}
        return A
    elseif A isa CuSparseMatrixCSC
        colptr = CuArray{Int32}(A.colPtr)
        rowval = CuArray{Int32}(A.rowVal)
        nzval = CuArray{ElType}(A.nzVal)
        return CuSparseMatrixCSC{ElType,Int32}(colptr, rowval, nzval, size(A))
    else
        colptr = CuArray{Int32}(A.colptr)
        rowval = CuArray{Int32}(A.rowval)
        nzval = CuArray{ElType}(A.nzval)
        return CuSparseMatrixCSC{ElType,Int32}(colptr, rowval, nzval, size(A))
    end
end
_convert_to_gpu_matrix(A::MatrixOperator, ElType) = MatrixOperator(_convert_to_gpu_matrix(A.A, ElType))
_convert_to_gpu_matrix(A::ScaledOperator, ElType) = ScaledOperator(A.λ, _convert_to_gpu_matrix(A.L, ElType))
_convert_to_gpu_matrix(A::AddedOperator, ElType) = AddedOperator(map(op -> _convert_to_gpu_matrix(op, ElType), A.ops))

# change the type of `ADOs` to match the type of HEOMLS matrix
_HandleVectorType(M::Type{<:CuSparseMatrixCSC}, V::SparseVector) = CuArray{_CType(eltype(M))}(V)

_HandleTraceVectorType(M::Type{<:CuSparseMatrixCSC}, V::SparseVector) = CuSparseVector{_CType(eltype(M))}(V)
end

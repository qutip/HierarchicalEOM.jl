module HierarchicalEOM_CUDAExt

using HierarchicalEOM
import HierarchicalEOM.HeomBase: AbstractHEOMLSMatrix, _Tr, _HandleVectorType
import QuantumToolbox: _CType
import CUDA
import CUDA: cu, CuArray
import CUDA.CUSPARSE: CuSparseVector, CuSparseMatrixCSC
import SparseArrays: sparse, SparseVector, SparseMatrixCSC

@doc raw"""
    cu(M::AbstractHEOMLSMatrix)
Return a new HEOMLS-matrix-type object with `M.data` is in the type of `CuSparseMatrixCSC{ComplexF32, Int32}` for gpu calculations.
"""
cu(M::AbstractHEOMLSMatrix) = CuSparseMatrixCSC(M)

@doc raw"""
    CuSparseMatrixCSC(M::AbstractHEOMLSMatrix)
Return a new HEOMLS-matrix-type object with `M.data` is in the type of `CuSparseMatrixCSC{ComplexF32, Int32}` for gpu calculations.
"""
function CuSparseMatrixCSC(M::T) where {T<:AbstractHEOMLSMatrix}
    A = M.data
    if typeof(A) <: CuSparseMatrixCSC
        return M
    else
        colptr = CuArray{Int32}(A.colptr)
        rowval = CuArray{Int32}(A.rowval)
        nzval = CuArray{ComplexF32}(A.nzval)
        A_gpu = CuSparseMatrixCSC{ComplexF32,Int32}(colptr, rowval, nzval, size(A))
        if T <: M_S
            return M_S(A_gpu, M.tier, M.dims, M.N, M.sup_dim, M.parity)
        elseif T <: M_Boson
            return M_Boson(A_gpu, M.tier, M.dims, M.N, M.sup_dim, M.parity, M.bath, M.hierarchy)
        elseif T <: M_Fermion
            return M_Fermion(A_gpu, M.tier, M.dims, M.N, M.sup_dim, M.parity, M.bath, M.hierarchy)
        else
            return M_Boson_Fermion(
                A_gpu,
                M.Btier,
                M.Ftier,
                M.dims,
                M.N,
                M.sup_dim,
                M.parity,
                M.Bbath,
                M.Fbath,
                M.hierarchy,
            )
        end
    end
end
function CuSparseMatrixCSC{ComplexF32}(M::HEOMSuperOp)
    A = M.data
    AType = typeof(A)
    if AType == CuSparseMatrixCSC{ComplexF32,Int32}
        return M
    elseif AType <: CuSparseMatrixCSC
        colptr = CuArray{Int32}(A.colPtr)
        rowval = CuArray{Int32}(A.rowVal)
        nzval = CuArray{ComplexF32}(A.nzVal)
        A_gpu = CuSparseMatrixCSC{ComplexF32,Int32}(colptr, rowval, nzval, size(A))
        return HEOMSuperOp(A_gpu, M.dims, M.N, M.parity)
    else
        colptr = CuArray{Int32}(A.colptr)
        rowval = CuArray{Int32}(A.rowval)
        nzval = CuArray{ComplexF32}(A.nzval)
        A_gpu = CuSparseMatrixCSC{ComplexF32,Int32}(colptr, rowval, nzval, size(A))
        return HEOMSuperOp(A_gpu, M.dims, M.N, M.parity)
    end
end

function _Tr(M::AbstractHEOMLSMatrix{T}) where {T<:CuSparseMatrixCSC}
    D = prod(M.dims)
    return CuSparseVector(SparseVector(M.N * D^2, [1 + n * (D + 1) for n in 0:(D-1)], ones(eltype(M), D)))
end

# change the type of `ADOs` to match the type of HEOMLS matrix
_HandleVectorType(::AbstractHEOMLSMatrix{<:CuSparseMatrixCSC{T}}, V::SparseVector) where {T<:Number} =
    CuArray{_CType(T)}(V)

##### We first remove this part because there are errors when solveing steady states using GPU
# function _HandleSteadyStateMatrix(MatrixType::Type{TM}, M::AbstractHEOMLSMatrix, S::Int) where TM <: CuSparseMatrixCSC
#     colptr = Vector{Int32}(M.data.colPtr)
#     rowval = Vector{Int32}(M.data.rowVal)
#     nzval  = Vector{ComplexF32}(M.data.nzVal)
#     A = SparseMatrixCSC{ComplexF32, Int32}(S, S, colptr, rowval, nzval)
#     A[1,1:S] .= 0f0
#     
#     # sparse(row_idx, col_idx, values, row_dims, col_dims)
#     A += sparse(ones(Int32, M.dim), [Int32((n - 1) * (M.dim + 1) + 1) for n in 1:(M.dim)], ones(ComplexF32, M.dim), S, S)
#     return CuSparseMatrixCSC(A)
# end

end

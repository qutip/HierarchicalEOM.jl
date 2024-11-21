module HierarchicalEOM_CUDAExt

using HierarchicalEOM
import HierarchicalEOM: _Tr, _HandleVectorType
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

_Tr(M::Type{<:CuSparseMatrixCSC}, dims::SVector, N::Int) = CuSparseVector(_Tr(eltype(M), dims, N))

# change the type of `ADOs` to match the type of HEOMLS matrix
_HandleVectorType(M::Type{<:CuSparseMatrixCSC}, V::SparseVector) = CuArray{_CType(eltype(M))}(V)

end

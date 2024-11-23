module HierarchicalEOM_CUDAExt

using HierarchicalEOM
import HierarchicalEOM: _Tr, _HandleVectorType
import QuantumToolbox: _CType
import CUDA
import CUDA: cu, CuArray
import CUDA.CUSPARSE: CuSparseVector, CuSparseMatrixCSC
import SparseArrays: sparse, SparseVector, SparseMatrixCSC
import SciMLOperators: MatrixOperator, ScaledOperator, AddedOperator

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
    A_gpu = _convert_to_gpu_matrix(M.data)
    if T <: M_S
        return M_S(A_gpu, M.tier, M.dims, M.N, M.sup_dim, M.parity)
    elseif T <: M_Boson
        return M_Boson(A_gpu, M.tier, M.dims, M.N, M.sup_dim, M.parity, M.bath, M.hierarchy)
    elseif T <: M_Fermion
        return M_Fermion(A_gpu, M.tier, M.dims, M.N, M.sup_dim, M.parity, M.bath, M.hierarchy)
    else
        return M_Boson_Fermion(A_gpu, M.Btier, M.Ftier, M.dims, M.N, M.sup_dim, M.parity, M.Bbath, M.Fbath, M.hierarchy)
    end
end

CuSparseMatrixCSC{ComplexF32}(M::HEOMSuperOp) = HEOMSuperOp(_convert_to_gpu_matrix(M.data), M.dims, M.N, M.parity)

function _convert_to_gpu_matrix(A::AbstractSparseMatrix)
    if A isa CuSparseMatrixCSC{ComplexF32,Int32}
        return A
    elseif A isa CuSparseMatrixCSC
        colptr = CuArray{Int32}(A.colPtr)
        rowval = CuArray{Int32}(A.rowVal)
        nzval = CuArray{ComplexF32}(A.nzVal)
        return CuSparseMatrixCSC{ComplexF32,Int32}(colptr, rowval, nzval, size(A))
    else
        colptr = CuArray{Int32}(A.colptr)
        rowval = CuArray{Int32}(A.rowval)
        nzval = CuArray{ComplexF32}(A.nzval)
        return CuSparseMatrixCSC{ComplexF32,Int32}(colptr, rowval, nzval, size(A))
    end
end
_convert_to_gpu_matrix(A::MatrixOperator) = MatrixOperator(_convert_to_gpu_matrix(A.A))
_convert_to_gpu_matrix(A::ScaledOperator) = ScaledOperator(A.λ, _convert_to_gpu_matrix(A.L))
_convert_to_gpu_matrix(A::AddedOperator) = AddedOperator(map(op -> _convert_to_gpu_matrix(op), A.ops))

_Tr(M::Type{<:CuSparseMatrixCSC}, dims::SVector, N::Int) = CuSparseVector(_Tr(eltype(M), dims, N))

# change the type of `ADOs` to match the type of HEOMLS matrix
_HandleVectorType(M::Type{<:CuSparseMatrixCSC}, V::SparseVector) = CuArray{_CType(eltype(M))}(V)

end

module HierarchicalEOM_CUDAExt

using HierarchicalEOM
import HierarchicalEOM.HeomAPI:  _HandleVectorType, _HandleSteadyStateMatrix
import HierarchicalEOM.Spectrum: _HandleIdentityType
import CUDA: cu, CuArray
import CUDA.CUSPARSE: CuSparseMatrixCSC
import SparseArrays: SparseVector

@doc raw"""
    cu(M::AbstractHEOMLSMatrix)
Return a new HEOMLS-matrix-type object with `M.data` is in the type of `CuSparseMatrixCSC{ComplexF32, Int32}` for gpu calculations.
"""
cu(M::AbstractHEOMLSMatrix) = CuSparseMatrixCSC(M)

@doc raw"""
    CuSparseMatrixCSC(M::AbstractHEOMLSMatrix)
Return a new HEOMLS-matrix-type object with `M.data` is in the type of `CuSparseMatrixCSC{ComplexF32, Int32}` for gpu calculations.
"""
function CuSparseMatrixCSC(M::T) where T <: AbstractHEOMLSMatrix
    A = M.data
    if typeof(A) <: CuSparseMatrixCSC
        return M
    else
        colptr = CuArray{Int32}(A.colptr)
        rowval = CuArray{Int32}(A.rowval)
        nzval  = CuArray{ComplexF32}(A.nzval)
        A_gpu  = CuSparseMatrixCSC{ComplexF32, Int32}(colptr, rowval, nzval, size(A))
        if T <: M_S
            return M_S(A_gpu, M.tier, M.dim, M.N, M.sup_dim, M.parity)
        elseif T <: M_Boson
            return M_Boson(A_gpu, M.tier, M.dim, M.N, M.sup_dim, M.parity, M.bath, M.hierarchy)
        elseif T <: M_Fermion
            return M_Fermion(A_gpu, M.tier, M.dim, M.N, M.sup_dim, M.parity, M.bath, M.hierarchy)
        else
            return M_Boson_Fermion(A_gpu, M.Btier, M.Ftier, M.dim, M.N, M.sup_dim, M.parity, M.Bbath, M.Fbath, M.hierarchy)
        end
    end
end

# for changing a `CuArray` back to `ADOs`
function _HandleVectorType(V::T, cp::Bool=false) where T <: CuArray
    return Vector{ComplexF64}(V)
end

# for changing the type of `ADOs` to match the type of HEOMLS matrix 
function _HandleVectorType(MatrixType::Type{TM}, V::SparseVector) where TM <: CuSparseMatrixCSC
    TE = eltype(MatrixType)
    return CuArray{TE}(V)
end

function _HandleSteadyStateMatrix(MatrixType::Type{TM}, M::AbstractHEOMLSMatrix, S::Int) where TM <: CuSparseMatrixCSC
    colptr = Vector{Int32}(M.data.colptr)
    rowval = Vector{Int32}(M.data.rowval)
    nzval  = Vector{ComplexF32}(M.data.rowval)
    A = SparseMatrixCSC{ComplexF32, Int32}(S, S, colptr, rowval, nzval)
    A[1,1:S] .= 0
    
    # sparse(row_idx, col_idx, values, row_dims, col_dims)
    A += sparse(ones(Int32, M.dim), [Int32((n - 1) * (M.dim + 1) + 1) for n in 1:(M.dim)], ones(ComplexF32, M.dim), S, S)
    return CuSparseMatrixCSC(A)
end

function _HandleIdentityType(MatrixType::Type{TM}, S::Int) where TM <: CuSparseMatrixCSC
    colptr = CuArray{Int32}(Int32(1):Int32(S+1))
    rowval = CuArray{Int32}(Int32(1):Int32(S))
    nzval  = CuArray{ComplexF32}(ones(ComplexF32, S))
    return CuSparseMatrixCSC{ComplexF32, Int32}(colptr, rowval, nzval, (S, S))
end

end
module HierarchicalEOM_CUDAExt

using HierarchicalEOM
import HierarchicalEOM:
    _reset_HEOMLS_data,
    _HandleVectorType,
    _HandleTraceVectorType,
    _HandleSteadyStateMatrix,
    _SteadyStateConstraint,
    _get_SciML_matrix_wrapper,
    get_cached_HEOMLS_data
import QuantumToolbox: _complex_float_type, _convert_eltype_wordsize, makeVal, getVal, get_typename_wrapper
import CUDA
import CUDA: cu, CuArray
import CUDA.CUSPARSE: CuSparseVector, CuSparseMatrixCSC, CuSparseMatrixCSR, AbstractCuSparseMatrix
import SparseArrays: AbstractSparseMatrix, sparse, SparseVector, SparseMatrixCSC
import LinearAlgebra: Diagonal
import SciMLOperators: MatrixOperator, ScaledOperator, AddedOperator, IdentityOperator, TensorProductOperator, AbstractSciMLOperator
import FillArrays: Eye

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

If `{T}` is not specified, default to `{eltype(M)}`.
"""
CuSparseMatrixCSC(M::MT) where {MT<:AbstractHEOMLSMatrix} =
    _reset_HEOMLS_data(M, _convert_to_gpu_matrix(M.data, CuSparseMatrixCSC{eltype(M)}))
CuSparseMatrixCSC{T}(M::MT) where {T,MT<:AbstractHEOMLSMatrix} =
    _reset_HEOMLS_data(M, _convert_to_gpu_matrix(M.data, CuSparseMatrixCSC{T}))

@doc raw"""
    CuSparseMatrixCSR{T}(M::AbstractHEOMLSMatrix)

Return a new HEOMLS-matrix-type object with `M.data` is in the type of `CUDA.CUSPARSE.CuSparseMatrixCSR` with element type `T` for gpu calculations.

If `{T}` is not specified, default to `{eltype(M)}`.
"""
CuSparseMatrixCSR(M::MT) where {MT<:AbstractHEOMLSMatrix} =
    _reset_HEOMLS_data(M, _convert_to_gpu_matrix(M.data, CuSparseMatrixCSR{eltype(M)}))
CuSparseMatrixCSR{T}(M::MT) where {T,MT<:AbstractHEOMLSMatrix} =
    _reset_HEOMLS_data(M, _convert_to_gpu_matrix(M.data, CuSparseMatrixCSR{T}))

CuSparseMatrixCSC{T}(M::HEOMSuperOp) where {T} =
    HEOMSuperOp(_convert_to_gpu_matrix(M.data, CuSparseMatrixCSC{T}), M.dimensions, M.N, M.parity)
CuSparseMatrixCSR{T}(M::HEOMSuperOp) where {T} =
    HEOMSuperOp(_convert_to_gpu_matrix(M.data, CuSparseMatrixCSR{T}), M.dimensions, M.N, M.parity)

_convert_to_gpu_matrix(A::AbstractSparseMatrix, MType::Type{T}) where {T<:AbstractCuSparseMatrix} = MType(A)
_convert_to_gpu_matrix(A::AbstractMatrix, MType::Type{T}) where {T<:AbstractCuSparseMatrix} = MType(sparse(A))

_convert_to_gpu_matrix(A::MatrixOperator, MType) = MatrixOperator(_convert_to_gpu_matrix(A.A, MType))
_convert_to_gpu_matrix(A::ScaledOperator, MType) = ScaledOperator(A.λ, _convert_to_gpu_matrix(A.L, MType))
_convert_to_gpu_matrix(A::AddedOperator, MType) = AddedOperator(map(op -> _convert_to_gpu_matrix(op, MType), A.ops))
_convert_to_gpu_matrix(A::IdentityOperator, MType) = A
_convert_to_gpu_matrix(A::TensorProductOperator, MType) =
    TensorProductOperator(_convert_to_gpu_matrix(A.ops[1], MType), _convert_to_gpu_matrix(A.ops[2], MType))

# change the type of `ADOs` to match the type of HEOMLS matrix
_HandleVectorType(M::Type{<:AbstractCuSparseMatrix}, V::SparseVector) = CuArray{_complex_float_type(eltype(M))}(V)

_HandleTraceVectorType(M::Type{<:AbstractCuSparseMatrix}, V::SparseVector) =
    CuSparseVector{_complex_float_type(eltype(M))}(V)

_HandleSteadyStateMatrix(
    M::AbstractHEOMLSMatrix{<:MatrixOperator{T,MT}},
    ::CuArray{T},
) where {T<:Number,MT<:AbstractCuSparseMatrix} =
    M.data.A + get_typename_wrapper(M.data.A)(_SteadyStateConstraint(T, prod(M.dimensions), size(M, 1)))

# if the HEOMLS Matrix is transferred to GPU by proper apis, it should have all the operators in the same sparse format
# To avoid scalar indexing in potential concretization, make the constraint the same type sparse format as M.data
# Do not specify the element type for the CuSparseMatrix... (No method)
_HandleSteadyStateMatrix(M::AbstractHEOMLSMatrix{<:AddedOperator{T}}, b::CuArray{T}) where {T<:Number} =
    get_cached_HEOMLS_data(
        M.data + _get_SciML_matrix_wrapper(M)(_SteadyStateConstraint(T, prod(M.dimensions), size(M, 1))),
        b,
    )
end

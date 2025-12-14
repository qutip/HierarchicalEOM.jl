module HierarchicalEOM_CUDAExt

using HierarchicalEOM
import HierarchicalEOM: _HandleVectorType, _HandleTraceVectorType, _HandleSteadyStateMatrix, _SteadyStateConstraint
import QuantumToolbox: _complex_float_type, _convert_eltype_wordsize, makeVal, getVal, get_typename_wrapper
import CUDA
import CUDA: cu, CuArray
import CUDA.CUSPARSE: CuSparseVector, CuSparseMatrixCSC, CuSparseMatrixCSR, AbstractCuSparseArray
import SparseArrays: AbstractSparseMatrix, sparse, SparseVector, SparseMatrixCSC
import SciMLOperators: MatrixOperator, ScaledOperator, AddedOperator, AbstractSciMLOperator

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
    _gen_gpu_HEOMLS(M, _convert_to_gpu_matrix(M.data, CuSparseMatrixCSC{eltype(M)}))
CuSparseMatrixCSC{T}(M::MT) where {T,MT<:AbstractHEOMLSMatrix} =
    _gen_gpu_HEOMLS(M, _convert_to_gpu_matrix(M.data, CuSparseMatrixCSC{T}))

@doc raw"""
    CuSparseMatrixCSR{T}(M::AbstractHEOMLSMatrix)

Return a new HEOMLS-matrix-type object with `M.data` is in the type of `CUDA.CUSPARSE.CuSparseMatrixCSR` with element type `T` for gpu calculations.

If `{T}` is not specified, default to `{eltype(M)}`.
"""
CuSparseMatrixCSR(M::MT) where {MT<:AbstractHEOMLSMatrix} =
    _gen_gpu_HEOMLS(M, _convert_to_gpu_matrix(M.data, CuSparseMatrixCSR{eltype(M)}))
CuSparseMatrixCSR{T}(M::MT) where {T,MT<:AbstractHEOMLSMatrix} =
    _gen_gpu_HEOMLS(M, _convert_to_gpu_matrix(M.data, CuSparseMatrixCSR{T}))

function _gen_gpu_HEOMLS(M::MT, data::AbstractSciMLOperator) where {MT<:AbstractHEOMLSMatrix}
    if M isa M_S
        return M_S(data, M.tier, M.dimensions, M.N, M.sup_dim, M.parity)
    elseif M isa M_Boson
        return M_Boson(data, M.tier, M.dimensions, M.N, M.sup_dim, M.parity, M.bath, M.hierarchy)
    elseif M isa M_Fermion
        return M_Fermion(data, M.tier, M.dimensions, M.N, M.sup_dim, M.parity, M.bath, M.hierarchy)
    else
        return M_Boson_Fermion(
            data,
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
    HEOMSuperOp(_convert_to_gpu_matrix(M.data, CuSparseMatrixCSC{T}), M.dimensions, M.N, M.parity)
CuSparseMatrixCSR{T}(M::HEOMSuperOp) where {T} =
    HEOMSuperOp(_convert_to_gpu_matrix(M.data, CuSparseMatrixCSR{T}), M.dimensions, M.N, M.parity)

_convert_to_gpu_matrix(A::AbstractSparseMatrix, MType::Type{T}) where {T<:AbstractCuSparseArray} = MType(A)
_convert_to_gpu_matrix(A::MatrixOperator, MType) = MatrixOperator(_convert_to_gpu_matrix(A.A, MType))
_convert_to_gpu_matrix(A::ScaledOperator, MType) = ScaledOperator(A.λ, _convert_to_gpu_matrix(A.L, MType))
_convert_to_gpu_matrix(A::AddedOperator, MType) = AddedOperator(map(op -> _convert_to_gpu_matrix(op, MType), A.ops))
_convert_to_gpu_matrix(A::TensorProductOperator, MType) = 
    TensorProductOperator(_convert_to_gpu_matrix(A.ops[1], MType), _convert_to_gpu_matrix(A.ops[2], MType))

# change the type of `ADOs` to match the type of HEOMLS matrix
_HandleVectorType(M::Type{<:AbstractCuSparseArray}, V::SparseVector) = CuArray{_complex_float_type(eltype(M))}(V)

_HandleTraceVectorType(M::Type{<:AbstractCuSparseArray}, V::SparseVector) =
    CuSparseVector{_complex_float_type(eltype(M))}(V)

_HandleSteadyStateMatrix(M::AbstractHEOMLSMatrix{<:MatrixOperator{T,MT}}) where {T<:Number,MT<:AbstractCuSparseArray} =
    M.data.A + get_typename_wrapper(M.data.A){eltype(M)}(_SteadyStateConstraint(T, prod(M.dimensions), size(M, 1)))
end

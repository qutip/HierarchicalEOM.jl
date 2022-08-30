abstract type AbstractHEOMMatrix end
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
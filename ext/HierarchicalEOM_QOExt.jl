module HierarchicalEOM_QOExt

import QuantumOptics: AbstractOperator, Operator
using HierarchicalEOM
import HierarchicalEOM.HeomBase: HandleMatrixType

@doc raw"""
    QOoperator(data::AbstractMatrix, refOP::AbstractOperator)
Return the operator under the type of `QuantumOptics.AbstractOperator` from a given matrix and the basis in reference operator.

# Parameters
- `data` : The matrix-type operator.
- `refOP<:QuantumOptics.AbstractOperator` : the reference operator which contains the basis defined in `QuantumOptics`.
"""
QOoperator(data::AbstractMatrix, refOP::AbstractOperator) = Operator(refOP.basis_l, refOP.basis_r, data)

function HandleMatrixType(M::AbstractOperator, dim::Int=0, MatrixName::String="")
    if dim > 0
        if size(M.data) == (dim, dim)
            return copy(M.data)
        else
            error("The size of matrix $(MatrixName) should be: ($(dim), $(dim)).")
        end
    elseif dim == 0
        N1, N2 = size(M.data)
        if N1 == N2
            return copy(M.data)
        else
            error("The size of matrix $(MatrixName) should be squared matrix.")
        end
    end
end

end
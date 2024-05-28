module HierarchicalEOM_QuantumToolboxExt

import QuantumToolbox: Qobj, QuantumObject
using HierarchicalEOM
import HierarchicalEOM.HeomBase: HandleMatrixType

@doc raw"""
    Qobj(data::AbstractMatrix, refOP::QuantumObject)
Return the operator under the type of `QuantumToolbox.QuantumObject` from a given matrix and reference operator (`type`, `dims`).

# Parameters
- `data` : The matrix-type operator.
- `refOP<:QuantumToolbox.QuantumObject` : the reference operator from `QuantumToolbox`.
"""
Qobj(data::AbstractMatrix, refOP::QuantumObject) = QuantumObject(data, refOP.type, refOP.dims)

@doc raw"""
    QuantumObject(data::AbstractMatrix, refOP::QuantumObject)
Return the operator under the type of `QuantumToolbox.QuantumObject` from a given matrix and reference operator (`type`, `dims`).

# Parameters
- `data` : The matrix-type operator.
- `refOP<:QuantumToolbox.QuantumObject` : the reference operator from `QuantumToolbox`.
"""
QuantumObject(data::AbstractMatrix, refOP::QuantumObject) = QuantumObject(data, refOP.type, refOP.dims)

function HandleMatrixType(M::QuantumObject, dim::Int = 0, MatrixName::String = "")
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

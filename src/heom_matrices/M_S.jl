export M_S

@doc raw"""
    struct M_S <: AbstractHEOMLSMatrix
HEOM Liouvillian superoperator matrix with cutoff level of the hierarchy equals to `0`.  
This corresponds to the standard Schrodinger (Liouville-von Neumann) equation, namely
```math
M[\cdot]=-i \left[H_{sys}, \cdot \right]_-,
```
where ``[\cdot, \cdot]_-`` stands for commutator.

# Fields
- `data<:AbstractSciMLOperator` : the matrix of HEOM Liouvillian superoperator
- `tier` : the tier (cutoff level) for the hierarchy, which equals to `0` in this case
- `dimensions` : the dimension list of the coupling operator (should be equal to the system dimensions).
- `N` : the number of total ADOs, which equals to `1` (only the reduced density operator) in this case
- `sup_dim` : the dimension of system superoperator
- `parity` : the parity label of the operator which HEOMLS is acting on (usually `EVEN`, only set as `ODD` for calculating spectrum of fermionic system).

!!! note "`dims` property"
    For a given `M::M_S`, `M.dims` or `getproperty(M, :dims)` returns its `dimensions` in the type of integer-vector.
"""
struct M_S{T<:AbstractSciMLOperator} <: AbstractHEOMLSMatrix{T}
    data::T
    tier::Int
    dimensions::Dimensions
    N::Int
    sup_dim::Int
    parity::AbstractParity
end

@doc raw"""
    M_S(Hsys, parity=EVEN; verbose=true)
Generate HEOM Liouvillian superoperator matrix with cutoff level of the hierarchy equals to `0`.  
This corresponds to the standard Schrodinger (Liouville-von Neumann) equation, namely
```math
M[\cdot]=-i \left[H_{sys}, \cdot \right]_-,
```
where ``[\cdot, \cdot]_-`` stands for commutator.

# Parameters
- `Hsys` : The time-independent system Hamiltonian or Liouvillian
- `parity::AbstractParity` : the parity label of the operator which HEOMLS is acting on (usually `EVEN`, only set as `ODD` for calculating spectrum of fermionic system).
- `verbose::Bool` : To display verbose output during the process or not. Defaults to `true`.

Note that the parity only need to be set as `ODD` when the system contains fermionic systems and you need to calculate the spectrum (density of states) of it.
"""
@noinline function M_S(Hsys::QuantumObject, parity::AbstractParity = EVEN; verbose::Bool = true)

    # check for system dimension
    _Hsys = HandleMatrixType(Hsys, "Hsys (system Hamiltonian or Liouvillian)")
    sup_dim = prod(_Hsys.dimensions)^2

    # the Liouvillian operator for free Hamiltonian
    if verbose
        println("Constructing Liouville-von Neumann superoperator...")
        flush(stdout)
    end
    Lsys = MatrixOperator(minus_i_L_op(_Hsys))
    if verbose
        println("[DONE]")
        flush(stdout)
    end
    return M_S(Lsys, 0, copy(_Hsys.dimensions), 1, sup_dim, parity)
end

_getBtier(M::M_S) = 0
_getFtier(M::M_S) = 0

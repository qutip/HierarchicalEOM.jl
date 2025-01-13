# [HEOMLS Matrix for Schrödinger Equation](@id doc-M_S)
The HEOM Liouvillian superoperator matrix with cutoff level of the hierarchy equals to `0`: [`struct M_S <: AbstractHEOMLSMatrix`](@ref M_S) 

This corresponds to the standard Schrodinger (Liouville-von Neumann) equation, namely
```math
\hat{\mathcal{M}}[\cdot]=-i \left[H_{s}, \cdot \right]_-,
```
where ``[\cdot, \cdot]_-`` stands for commutator.

## Construct Matrix
To construct the HEOM matrix for Schrödinger Equation, one can call 

[`M_S(Hsys, parity)`](@ref M_S) with the following parameters:

*args* (Arguments)
 - `Hsys` : The time-independent system Hamiltonian
 - `parity::AbstractParity` : the [parity](@ref doc-Parity) label of the operator which HEOMLS is acting on. Defaults to `EVEN`.
*kwargs* (Keyword Arguments)
 - `verbose::Bool` : To display verbose output during the process or not. Defaults to `true`.

For example:
```julia
Hs::QuantumObject # system Hamiltonian

# create HEOMLS matrix in both EVEN and ODD parity
M_even = M_S(Hs) 
M_odd  = M_S(Hs, ODD) 
```

## Fields
The fields of the structure [`M_S`](@ref) are as follows:
 - `data` : the sparse matrix of HEOM Liouvillian superoperator
 - `tier` : the tier (cutoff level) for the hierarchy, which equals to `0` in this case
 - `dimensions` : the dimension list of the coupling operator (should be equal to the system dimensions).
 - `N` : the number of total [ADOs](@ref doc-ADOs), which equals to `1` (only the reduced density operator) in this case
 - `sup_dim` : the dimension of system superoperator
 - `parity::AbstractParity` : the [parity](@ref doc-Parity) label of the operator which HEOMLS is acting on.

One can obtain the value of each fields as follows:
```julia
M::M_S

M.data
M.tier
M.dimensions
M.dims
M.N
M.sup_dim
M.parity
```
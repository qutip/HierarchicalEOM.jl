# [HEOMLS Matrix for Fermionic Baths](@id doc-M_Fermion)
The HEOM Liouvillian superoperator matrix [`struct M_Fermion <: AbstractHEOMLSMatrix`](@ref M_Fermion) which describes the interactions between the system and multiple [Fermionic baths](@ref doc-Fermionic-Bath).

## Construct Matrix
To construct the HEOM matrix in this case, one can call 

[`M_Fermion(Hsys, tier, Bath, parity)`](@ref M_Fermion) with the following parameters:

*args* (Arguments)
 - `Hsys` : The time-independent system Hamiltonian
 - `tier::Int` : the tier (cutoff level) for the fermionic bath
 - `Bath::Vector{FermionBath}` : objects for different [fermionic baths](@ref doc-Fermionic-Bath)
 - `parity::AbstractParity` : the [parity](@ref doc-Parity) label of the operator which HEOMLS is acting on. Defaults to `EVEN`.

*kwargs* (Keyword Arguments)
 - `threshold::Real` : The threshold of the [importance value](@ref doc-Importance-Value-and-Threshold). Defaults to `0.0`.
 - `verbose::Bool` : To display verbose output and progress bar during the process or not. Defaults to `true`.

For example:
```julia
Hs::QuantumObject # system Hamiltonian
tier = 3
Bath::FermionBath

# create HEOMLS matrix in both EVEN and ODD parity
M_even = M_Fermion(Hs, tier, Bath) 
M_odd  = M_Fermion(Hs, tier, Bath, ODD) 
```

## Fields
The fields of the structure [`M_Fermion`](@ref) are as follows:
 - `data` : the sparse matrix of HEOM Liouvillian superoperator
 - `tier` : the tier (cutoff level) for the fermionic hierarchy
 - `dimensions` : the dimension list of the coupling operator (should be equal to the system dimensions).
 - `N` : the number of total [ADOs](@ref doc-ADOs)
 - `sup_dim` : the dimension of system superoperator
 - `parity` : the [parity](@ref doc-Parity) label of the operator which HEOMLS is acting on.
 - `bath::Vector{FermionBath}` : the vector which stores all [`FermionBath`](@ref doc-Fermionic-Bath) objects
 - `hierarchy::HierarchyDict`: the object which contains all [dictionaries](@ref doc-Hierarchy-Dictionary) for fermion-bath-ADOs hierarchy.

One can obtain the value of each fields as follows:
```julia
M::M_Fermion

M.data
M.tier
M.dimensions
M.dims
M.N
M.sup_dim
M.parity
M.bath
M.hierarchy
```
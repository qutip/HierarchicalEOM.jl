# [HEOMLS Matrix for Fermionic Baths](@id doc-M_Fermion)
The HEOM Liouvillian superoperator matrix [`struct M_Fermion <: AbstractHEOMMatrix`](@ref M_Fermion) which describes the interactions between the system and multiple [Fermionic baths](@ref doc-Fermionic-Bath).

## Construct Matrix
To construct the HEOM matrix in this case, one can call 

[`M_Fermion(Hsys, tier, Bath, parity)`](@ref M_Fermion_Constructor) with the following parameters:

*args* (Arguments)
 - `Hsys` : The time-independent system Hamiltonian
 - `tier::Int` : the tier (cutoff level) for the fermionic bath
 - `Bath::Vector{FermionBath}` : objects for different [fermionic baths](@ref doc-Fermionic-Bath)
 - `parity::Symbol` : the [parity](@ref doc-Parity) label of the fermionic system. Defaults to `:even`.

*kwargs* (Keyword Arguments)
 - `threshold::Real` : The threshold of the [importance value](@ref doc-Importance-Value-and-Threshold). Defaults to `0.0`.
 - `verbose::Bool` : To display verbose output and progress bar during the process or not. Defaults to `true`.

For example:
```julia
Hs::AbstractMatrix # system Hamiltonian
tier = 3
Bath::FermionBath

# create HEOMLS matrix in both :even and :odd parity
M_even = M_Fermion(Hs, tier, Bath) 
M_odd  = M_Fermion(Hs, tier, Bath, :odd) 
```

## Fields
The fields of the structure [`M_Fermion`](@ref) are as follows:
 - `data` : the sparse matrix of HEOM Liouvillian superoperator
 - `tier` : the tier (cutoff level) for the fermionic hierarchy
 - `dim` : the dimension of system
 - `N` : the number of total [ADOs](@ref doc-ADOs)
 - `sup_dim` : the dimension of system superoperator
 - `parity` : the [parity](@ref doc-Parity) label of the fermionic system.
 - `bath::Vector{FermionBath}` : the vector which stores all [`FermionBath`](@ref doc-Fermionic-Bath) objects
 - `hierarchy::HierarchyDict`: the object which contains all [dictionaries](@ref doc-Hierarchy-Dictionary) for fermion-bath-ADOs hierarchy.

One obtain the value of each fields as follows:
```julia
M::M_Fermion

M.data
M.tier
M.dim
M.N
M.sup_dim
M.parity
M.bath
M.hierarchy
```
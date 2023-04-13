# [HEOMLS Matrix for Hybrid (Bosonic and Fermionic) Baths](@id doc-M_Boson_Fermion)
The HEOM Liouvillian superoperator matrix [`struct M_Boson_Fermion <: AbstractHEOMMatrix`](@ref M_Boson_Fermion) which describes the system simultaneously interacts with multiple [Bosonic baths](@ref doc-Bosonic-Bath) and [Fermionic baths](@ref doc-Fermionic-Bath). 

## Construct Matrix
To construct the HEOM matrix in this case, one can call 

[`M_Boson_Fermion(Hsys, Btier, Ftier, Bbath, Fbath, parity)`](@ref M_Boson_Fermion_Constructor) with the following parameters:

*args* (Arguments)
 - `Hsys` : The time-independent system Hamiltonian
 - `Btier::Int` : the tier (cutoff level) for the bosonic bath
 - `Ftier::Int` : the tier (cutoff level) for the fermionic bath
 - `Bbath::Vector{BosonBath}` : objects for different [bosonic baths](@ref doc-Bosonic-Bath)
 - `Fbath::Vector{FermionBath}` : objects for different [fermionic baths](@ref doc-Fermionic-Bath)
 - `parity::Symbol` : the [parity](@ref doc-Parity) label of the fermionic system. Defaults to `:even`.

*kwargs* (Keyword Arguments)
 - `threshold::Real` : The threshold of the [importance value](@ref doc-Importance-Value-and-Threshold). Defaults to `0.0`.
 - `verbose::Bool` : To display verbose output and progress bar during the process or not. Defaults to `true`.

For example:
```julia
Hs::AbstractMatrix # system Hamiltonian
Btier = 3
Ftier = 4
Bbath::BosonBath
Fbath::FermionBath

# create HEOMLS matrix in both :even and :odd parity
M_even = M_Fermion(Hs, Btier, Ftier, Bbath, Fbath) 
M_odd  = M_Fermion(Hs, Btier, Ftier, Bbath, Fbath, :odd) 
```

## Fields
The fields of the structure [`M_Boson_Fermion`](@ref) are as follows:
 - `data` : the sparse matrix of HEOM liouvillian superoperator
 - `Btier` : the tier (cutoff level) for bosonic hierarchy
 - `Ftier` : the tier (cutoff level) for fermionic hierarchy
 - `dim` : the dimension of system
 - `N` : the number of total [ADOs](@ref doc-ADOs)
 - `sup_dim` : the dimension of system superoperator
 - `parity` : the [parity](@ref doc-Parity) label of the fermionic system.
 - `Bbath::Vector{BosonBath}` : the vector which stores all [`BosonBath`](@ref doc-Bosonic-Bath) objects
 - `Fbath::Vector{FermionBath}` : the vector which stores all [`FermionBath`](@ref doc-Fermionic-Bath) objects
 - `hierarchy::MixHierarchyDict`: the object which contains all [dictionaries](@ref doc-Hierarchy-Dictionary) for mixed-bath-ADOs hierarchy.

One obtain the value of each fields as follows:
```julia
M::M_Boson_Fermion

M.data
M.Btier
M.Ftier
M.dim
M.N
M.sup_dim
M.parity
M.Bbath
M.Fbath
M.hierarchy
```
# [HEOMLS Matrix for Bosonic Baths](@id doc-M_Boson)
The HEOM Liouvillian superoperator matrix [`struct M_Boson <: AbstractHEOMLSMatrix`](@ref M_Boson) which describes the interactions between the system and multiple [Bosonic baths](@ref doc-Bosonic-Bath).

## Construct Matrix
To construct the HEOM matrix in this case, one can call 

[`M_Boson(Hsys, tier, Bath, parity)`](@ref M_Boson) with the following parameters:

*args* (Arguments)
 - `Hsys` : The time-independent system Hamiltonian
 - `tier::Int` : the tier (cutoff level) for the bosonic bath
 - `Bath::Vector{BosonBath}` : objects for different [bosonic baths](@ref doc-Bosonic-Bath)
 - `parity::AbstractParity` : the [parity](@ref doc-Parity) label of the operator which HEOMLS is acting on. Defaults to `EVEN`.

*kwargs* (Keyword Arguments)
 - `threshold::Real` : The threshold of the [importance value](@ref doc-Importance-Value-and-Threshold). Defaults to `0.0`.
 - `verbose::Bool` : To display verbose output and progress bar during the process or not. Defaults to `true`.

For example:
```julia
Hs::AbstractMatrix # system Hamiltonian
tier = 3
Bath::BosonBath

# create HEOMLS matrix in both :even and :odd parity
M_even = M_Boson(Hs, tier, Bath) 
M_odd  = M_Boson(Hs, tier, Bath, :odd) 
```

## Fields
The fields of the structure [`M_Boson`](@ref) are as follows:
 - `data` : the sparse matrix of HEOM Liouvillian superoperator
 - `tier` : the tier (cutoff level) for the bosonic hierarchy
 - `dim` : the dimension of system
 - `N` : the number of total [ADOs](@ref doc-ADOs)
 - `sup_dim` : the dimension of system superoperator
 - `parity` : the [parity](@ref doc-Parity) label of the operator which HEOMLS is acting on.
 - `bath::Vector{BosonBath}` : the vector which stores all [`BosonBath`](@ref doc-Bosonic-Bath) objects
 - `hierarchy::HierarchyDict`: the object which contains all [dictionaries](@ref doc-Hierarchy-Dictionary) for boson-bath-ADOs hierarchy.

One obtain the value of each fields as follows:
```julia
M::M_Boson

M.data
M.tier
M.dim
M.N
M.sup_dim
M.parity
M.bath
M.hierarchy
```
# HEOM API Module

## HEOM Liouvillian superoperator matrices

### HEOM matrix for standard Schrodinger (Liouville-von Neumann) equation 
#### Type Definition
```@docs
M_S
```
#### [Constructor](@id M_S_Constructor)
```@docs
M_S(Hsys, parity::Symbol=:even; verbose::Bool=true)
```

### Boson HEOM matrix
#### Type Definition
```@docs
M_Boson
```
#### [Constructor](@id M_Boson_Constructor)
```@docs
M_Boson(Hsys, tier::Int, Bath::Vector{BosonBath}, parity::Symbol=:even; threshold::Real=0.0, verbose::Bool=true)
```

### Fermion HEOM matrix
#### Type Definition
```@docs
M_Fermion
```
#### [Constructor](@id M_Fermion_Constructor)
```@docs
M_Fermion(Hsys, tier::Int, Bath::Vector{FermionBath}, parity::Symbol=:even; threshold::Real=0.0, verbose::Bool=true)
```

### Boson-Fermion HEOM matrix
#### Type Definition
```@docs
M_Boson_Fermion
```
#### [Constructor](@id M_Boson_Fermion_Constructor)
```@docs
M_Boson_Fermion(Hsys, tier_b::Int, tier_f::Int, Bath_b::Vector{BosonBath}, Bath_f::Vector{FermionBath}, parity::Symbol=:even; threshold::Real=0.0, verbose::Bool=true)
```

### Functions for HEOM Matrices
```@docs
size(M::AbstractHEOMMatrix)
Propagator
addBosonDissipator
addFermionDissipator
addTerminator
```

## Auxiliary Density Operators (ADOs)
```@docs
ADOs
```

```@docs
ADOs(V::AbstractVector, N::Int)
```

### Functions for Auxiliary Density Operators
```@docs
length(A::ADOs)
getRho
getADO(ados::ADOs, idx::Int)
Expect(op, ados::ADOs; take_real=true)
Expect(op, ados_list::Vector{ADOs}; take_real=true)
```

## [Hierarchy Dictionary](@id lib-Hierarchy-Dictionary)
```@docs
Nvec
```

```@docs
HierarchyDict
MixHierarchyDict
```

```@docs
getIndexEnsemble
```
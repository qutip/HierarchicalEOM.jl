# Heom API Module

## Heom liouvillian superoperator matrices

### Heom matrix for standard Schrodinger (Liouville-von Neumann) equation 
```@docs
M_S
```

```@docs
M_S(Hsys, parity::Symbol=:even; verbose::Bool=true)
```

### Boson Heom matrix
```@docs
M_Boson
```

```@docs
M_Boson(Hsys, tier::Int, Bath::Vector{BosonBath}, parity::Symbol=:even; threshold::Real=0.0, verbose::Bool=true)
```

### Fermion Heom matrix
```@docs
M_Fermion
```

```@docs
M_Fermion(Hsys, tier::Int, Bath::Vector{FermionBath}, parity::Symbol=:even; threshold::Real=0.0, verbose::Bool=true)
```

### Boson-Fermion Heom matrix
```@docs
M_Boson_Fermion
```

```@docs
M_Boson_Fermion(Hsys, tier_b::Int, tier_f::Int, Bath_b::Vector{BosonBath}, Bath_f::Vector{FermionBath}, parity::Symbol=:even; threshold::Real=0.0, verbose::Bool=true)
```

### Functions for Heom Matrices
```@docs
size(M::AbstractHEOMMatrix)
Propagator
addBosonicDissipator
addFermionicDissipator
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
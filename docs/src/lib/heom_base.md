# Heom Base Module API

## Heom liouvillian superoperator matrices

### Boson Heom matrix
```@docs
M_Boson
```

```@docs
M_Boson(Hsys::AbstractMatrix, tier::Int, Bath::Vector{BosonBath}; progressBar::Bool=true)
```

### Fermion Heom matrix
```@docs
M_Fermion
```

```@docs
M_Fermion(Hsys::AbstractMatrix, tier::Int, Bath::Vector{FermionBath}, parity::Symbol=:even; progressBar::Bool=true)
```

### Boson-Fermion Heom matrix
```@docs
M_Boson_Fermion
```

```@docs
M_Boson_Fermion(Hsys::AbstractMatrix, tier_b::Int, tier_f::Int, Bath_b::Vector{BosonBath}, Bath_f::Vector{FermionBath}, parity::Symbol=:even; progressBar::Bool=true)
```

### Functions for Heom Matrices
```@docs
size(M::AbstractHEOMMatrix)
addDissipator!
addTerminator!
```

## Auxiliary Density Operators (ADOs)
```@docs
ADOs
```

```@docs
ADOs(V::AbstractVector, dim::Int, Nb::Int=0, Nf::Int=0)
```

### Functions for Auxiliary Density Operators
```@docs
size(A::ADOs)
getRho
getADO(ados::ADOs, idx::Int)
getADO(ados::ADOs, idx_b::Int, idx_f::Int)
```
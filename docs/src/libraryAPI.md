# Library API

## Contents

```@contents
Pages = ["libraryAPI.md"]
```

## Index

```@index
Pages = ["libraryAPI.md"]
```

## [Bath Module](@id lib-Bath)

```@docs
C(bath::BosonBath, tlist::AbstractVector)
C(bath::FermionBath, tlist::AbstractVector)
Exponent
BosonBath
BosonBath(op::AbstractMatrix, η::Vector{Ti}, γ::Vector{Tj}, δ::Number=0.0; combine::Bool=true) where {Ti, Tj <: Number}
BosonBath(op::AbstractMatrix, η_real::Vector{Ti}, γ_real::Vector{Tj}, η_imag::Vector{Tk}, γ_imag::Vector{Tl}, δ::Tm=0.0; combine::Bool=true) where {Ti, Tj, Tk, Tl, Tm <: Number}
bosonReal
bosonReal(op::AbstractMatrix, η_real::Vector{Ti}, γ_real::Vector{Tj}) where {Ti, Tj <: Number}
bosonImag
bosonImag(op::AbstractMatrix, η_real::Vector{Ti}, γ_real::Vector{Tj}) where {Ti, Tj <: Number}
bosonRealImag
bosonRealImag(op::AbstractMatrix, η_real::Vector{Ti}, η_imag::Vector{Tj}, γ::Vector{Tk}) where {Ti, Tj, Tk <: Number}
BosonBathRWA
bosonAbsorb
bosonAbsorb(op, η_absorb::Vector{Ti}, γ_absorb::Vector{Tj}, η_emit::Vector{Tk}) where {Ti, Tj, Tk <: Number}
bosonEmit
bosonEmit(op, η_emit::Vector{Ti}, γ_emit::Vector{Tj}, η_absorb::Vector{Tk}) where {Ti, Tj, Tk <: Number}
FermionBath
FermionBath(op::AbstractMatrix, η_absorb::Vector{Ti}, γ_absorb::Vector{Tj}, η_emit::Vector{Tk}, γ_emit::Vector{Tl}, δ::Tm=0.0) where {Ti, Tj, Tk, Tl, Tm <: Number}
fermionAbsorb
fermionAbsorb(op::AbstractMatrix, η_absorb::Vector{Ti}, γ_absorb::Vector{Tj}, η_emit::Vector{Tk}) where {Ti, Tj, Tk <: Number}
fermionEmit
fermionEmit(op::AbstractMatrix, η_emit::Vector{Ti}, γ_emit::Vector{Tj}, η_absorb::Vector{Tk}) where {Ti, Tj, Tk <: Number}
```

## Correlation Functions

```@docs
Boson_DrudeLorentz_Matsubara
Boson_DrudeLorentz_Pade
Fermion_Lorentz_Matsubara
Fermion_Lorentz_Pade
```

## HEOM Liouvillian superoperator matrices
```@docs
M_S
M_S(Hsys, parity::Symbol=:even; verbose::Bool=true)
M_Boson
M_Boson(Hsys, tier::Int, Bath::Vector{BosonBath}, parity::Symbol=:even; threshold::Real=0.0, verbose::Bool=true)
M_Fermion
M_Fermion(Hsys, tier::Int, Bath::Vector{FermionBath}, parity::Symbol=:even; threshold::Real=0.0, verbose::Bool=true)
M_Boson_Fermion
M_Boson_Fermion(Hsys, tier_b::Int, tier_f::Int, Bath_b::Vector{BosonBath}, Bath_f::Vector{FermionBath}, parity::Symbol=:even; threshold::Real=0.0, verbose::Bool=true)
size(M::AbstractHEOMLSMatrix)
size(M::AbstractHEOMLSMatrix, dim::Int)
eltype(M::AbstractHEOMLSMatrix)
Propagator
addBosonDissipator
addFermionDissipator
addTerminator
```

## Auxiliary Density Operators (ADOs)
```@docs
ADOs
ADOs(V::AbstractVector, N::Int)
length(A::ADOs)
eltype(A::ADOs)
getRho
getADO
Expect
```

## [Hierarchy Dictionary](@id lib-Hierarchy-Dictionary)
```@docs
Nvec
HierarchyDict
MixHierarchyDict
getIndexEnsemble
```

## [Time Evolution](@id lib-Time-Evolution)
There are six function definitions of `evolution`, which depend on different input types and methods to solve the time evolution:
```@docs
evolution
```

## Steady State
There are three function definitions of `SteadyState`, which depend on different input types and methods to solve the stationary state:
```@docs
SteadyState
```

## Spectrum
```@docs
spectrum
```

## Misc.
```@docs
HierarchicalEOM.versioninfo()
```
The outputs will be something like the following:
```@example
using HierarchicalEOM
HierarchicalEOM.versioninfo()
```

```@docs
HierarchicalEOM.print_logo(io::IO=stdout)
```
The output will be something like the following:
```@example
using HierarchicalEOM
HierarchicalEOM.print_logo()
```
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
correlation_function
Exponent
BosonBath
BosonBath(op::QuantumObject, η::Vector{Ti}, γ::Vector{Tj}, δ::Number=0.0; combine::Bool=true) where {Ti <: Number, Tj <: Number}
BosonBath(op::QuantumObject, η_real::Vector{Ti}, γ_real::Vector{Tj}, η_imag::Vector{Tk}, γ_imag::Vector{Tl}, δ::Tm=0.0; combine::Bool=true) where {Ti <: Number, Tj <: Number, Tk <: Number, Tl <: Number, Tm <: Number}
bosonReal
bosonReal(op::QuantumObject, η_real::Vector{Ti}, γ_real::Vector{Tj}) where {Ti <: Number, Tj <: Number}
bosonImag
bosonImag(op::QuantumObject, η_real::Vector{Ti}, γ_real::Vector{Tj}) where {Ti <: Number, Tj <: Number}
bosonRealImag
bosonRealImag(op::QuantumObject, η_real::Vector{Ti}, η_imag::Vector{Tj}, γ::Vector{Tk}) where {Ti <: Number, Tj <: Number, Tk <: Number}
BosonBathRWA
bosonAbsorb
bosonAbsorb(op::QuantumObject, η_absorb::Vector{Ti}, γ_absorb::Vector{Tj}, η_emit::Vector{Tk}) where {Ti <: Number, Tj <: Number, Tk <: Number}
bosonEmit
bosonEmit(op::QuantumObject, η_emit::Vector{Ti}, γ_emit::Vector{Tj}, η_absorb::Vector{Tk}) where {Ti <: Number, Tj <: Number, Tk <: Number}
FermionBath
FermionBath(op::QuantumObject, η_absorb::Vector{Ti}, γ_absorb::Vector{Tj}, η_emit::Vector{Tk}, γ_emit::Vector{Tl}, δ::Tm=0.0) where {Ti <: Number, Tj <: Number, Tk <: Number, Tl <: Number, Tm <: Number}
fermionAbsorb
fermionAbsorb(op::QuantumObject, η_absorb::Vector{Ti}, γ_absorb::Vector{Tj}, η_emit::Vector{Tk}) where {Ti <: Number, Tj <: Number, Tk <: Number}
fermionEmit
fermionEmit(op::QuantumObject, η_emit::Vector{Ti}, γ_emit::Vector{Tj}, η_absorb::Vector{Tk}) where {Ti <: Number, Tj <: Number, Tk <: Number}
```

## Bath Correlation Functions

```@docs
Boson_DrudeLorentz_Matsubara
Boson_DrudeLorentz_Pade
Boson_Underdamped_Matsubara
Fermion_Lorentz_Matsubara
Fermion_Lorentz_Pade
```

## Parity
```@docs
EvenParity
OddParity
EVEN
ODD
```

## HEOM Liouvillian superoperator matrices
```@docs
HEOMSuperOp
HEOMSuperOp(op, opParity::AbstractParity, refHEOMLS::AbstractHEOMLSMatrix)
HEOMSuperOp(op, opParity::AbstractParity, refADOs::ADOs)
HEOMSuperOp(op, opParity::AbstractParity, dims, N::Int)
AbstractHEOMLSMatrix
M_S
M_S(Hsys::QuantumObject, parity::AbstractParity=EVEN; verbose::Bool=true)
M_Boson
M_Boson(Hsys::QuantumObject, tier::Int, Bath::Vector{BosonBath}, parity::AbstractParity=EVEN; threshold::Real=0.0, assemble::Union{Val,Symbol}=Val(:full), verbose::Bool=true)
M_Fermion
M_Fermion(Hsys::QuantumObject, tier::Int, Bath::Vector{FermionBath}, parity::AbstractParity=EVEN; threshold::Real=0.0, assemble::Union{Val,Symbol}=Val(:full), verbose::Bool=true)
M_Boson_Fermion
M_Boson_Fermion(Hsys::QuantumObject, tier_b::Int, tier_f::Int, Bath_b::Vector{BosonBath}, Bath_f::Vector{FermionBath}, parity::AbstractParity=EVEN; threshold::Real=0.0, assemble::Union{Val,Symbol}=Val(:full), verbose::Bool=true)
Base.size(M::HEOMSuperOp)
Base.size(M::HEOMSuperOp, dim::Int)
Base.size(M::AbstractHEOMLSMatrix)
Base.size(M::AbstractHEOMLSMatrix, dim::Int)
Base.eltype(M::HEOMSuperOp)
Base.eltype(M::AbstractHEOMLSMatrix)
iscached
isconstant
SciMLOperators.concretize(M::AbstractHEOMLSMatrix)
Propagator
addBosonDissipator
addFermionDissipator
addTerminator
```

## Auxiliary Density Operators (ADOs)
```@docs
ADOs
ADOs(V::AbstractVector, N::Int)
Base.length(A::ADOs)
Base.eltype(A::ADOs)
getRho
getADO
QuantumToolbox.expect
```

## [Hierarchy Dictionary](@id lib-Hierarchy-Dictionary)
```@docs
Nvec
HierarchyDict
MixHierarchyDict
getIndexEnsemble
```

## [Time Evolution](@id lib-Time-Evolution)
There are two function definitions of `HEOMsolve`, which depend on different methods to solve the time evolution:
```@docs
HEOMsolveProblem
HEOMsolve
HEOMsolve_map
TimeEvolutionHEOMSol
```

## Stationary State
There are two function definitions of `steadystate`, which depend on different methods to solve the stationary state:
```@docs
steadystate
```

## Correlation Functions and Spectrum
```@docs
correlation_3op_2t
correlation_3op_1t
correlation_2op_2t
correlation_2op_1t
PowerSpectrum
DensityOfStates
```

## Misc.
```@docs
HierarchicalEOM.versioninfo()
HierarchicalEOM.about()
HierarchicalEOM.cite()
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
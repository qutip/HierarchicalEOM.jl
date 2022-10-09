# Bath Module API

## Bosonic Bath
```@docs
BosonBath
```

```@docs
BosonBath(op::AbstractMatrix, η::Vector{Ti}, γ::Vector{Tj}, δ::Number=0.0; combine::Bool=true) where {Ti, Tj <: Number}
BosonBath(op::AbstractMatrix, η_real::Vector{Ti}, γ_real::Vector{Tj}, η_imag::Vector{Tk}, γ_imag::Vector{Tl}, δ::Tm=0.0; combine::Bool=true) where {Ti, Tj, Tk, Tl, Tm <: Number}
```

```@docs
bosonReal
```

```@docs
bosonReal(op::AbstractMatrix, η_real::Vector{Ti}, γ_real::Vector{Tj}) where {Ti, Tj <: Number}
```

```@docs
bosonImag
```

```@docs
bosonImag(op::AbstractMatrix, η_real::Vector{Ti}, γ_real::Vector{Tj}) where {Ti, Tj <: Number}
```

```@docs
bosonRealImag
```

```@docs
bosonRealImag(op::AbstractMatrix, η_real::Vector{Ti}, η_imag::Vector{Tj}, γ::Vector{Tk}) where {Ti, Tj, Tk <: Number}
```

## Fermionic Bath
```@docs
FermionBath
```

```@docs
FermionBath(op::AbstractMatrix, η_absorb::Vector{Ti}, γ_absorb::Vector{Tj}, η_emit::Vector{Tk}, γ_emit::Vector{Tl}, δ::Tm=0.0) where {Ti, Tj, Tk, Tl, Tm <: Number}
```

```@docs
fermionAbsorb
```

```@docs
fermionAbsorb(op::AbstractMatrix, η_absorb::Vector{Ti}, γ_absorb::Vector{Tj}, η_emit::Vector{Tk}) where {Ti, Tj, Tk <: Number}
```

```@docs
fermionEmit
```

```@docs
fermionEmit(op::AbstractMatrix, η_emit::Vector{Ti}, γ_emit::Vector{Tj}, η_absorb::Vector{Tk}) where {Ti, Tj, Tk <: Number}
```

## Exponent
```@docs
Exponent
```
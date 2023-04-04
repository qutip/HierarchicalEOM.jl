# [Spectrum](@id doc-Spectrum)
We briefly summarize how to numerically compute the spectrum associated with the system degree of freedom. [Phys. Rev. Lett. 109, 266403 (2012)](https://link.aps.org/doi/10.1103/PhysRevLett.109.266403) showed that the spectrum can be evaluated either in time or frequency domain.

`Heom.jl` provides a function [`spectrum`](@ref) which performs the calculation in frequency domain. There are two different methods (as shown below) which depends on the [parity](@ref doc-Parity) of the HEOMLS matrices ``\hat{\mathcal{M}}`` corresponds to different system degree of freedom. 

If you want to calculate the spectrum associated with
 - **bosonic systems (Power Spectral Density)**: you have to provide ``\hat{\mathcal{M}}`` constructed in `:even` parity.
 - **fermionic systems (Density of States)**: you have to provide ``\hat{\mathcal{M}}`` constructed in `:odd` parity.

The function [`spectrum`](@ref) will automatically detect the [parity](@ref doc-Parity) of ``\hat{\mathcal{M}}`` by itself. Furthermore, the output of the function [`spectrum`](@ref) for both cases will always be in the type of `Vector{Float64}`, which contains the list of the spectrum values corresponding to the given `ω_list`.

See also the docstring : [`spectrum`](@ref)

## Power Spectral Density

```julia
M::AbstractHEOMMatrix # need to be in ":even" parity

# the input state can be in either type (but usually ADOs):
ρ::AbstractMatrix # the reduced density operator
ρ::ADOs # the ADOs solved from "evolution" or "SteadyState"

# the (usually annihilation) operator "a" as shown above
a::AbstractMatrix 

# the spectrum value for the specific frequency ω which need to be solved
ω_list = 0:0.5:2 # [0.0, 0.5, 1.0, 1.5, 2.0]

πSω = spectrum(M, ρ, a, ω_list)
```

## Density of States

```julia
M::AbstractHEOMMatrix # need to be in ":odd" parity

# the input state can be in either type (but usually ADOs):
ρ::AbstractMatrix # the reduced density operator
ρ::ADOs # the ADOs solved from "evolution" or "SteadyState"

# the (usually annihilation) operator "d" as shown above
d::AbstractMatrix 

# the spectrum value for the specific frequency ω which need to be solved
ω_list = 0:0.5:2 # [0.0, 0.5, 1.0, 1.5, 2.0]

πAω = spectrum(M, ρ, d, ω_list)
```
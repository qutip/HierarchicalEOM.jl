# Physical Analysis Functions

## Time Evolution
Using the method based on generating [`Propagator`](@ref)
```@docs
evolution(M::AbstractHEOMMatrix, ρ0, Δt::Real, steps::Int; threshold   = 1.0e-6, nonzero_tol = 1.0e-14, verbose::Bool = true, filename::String = "")
evolution(M::AbstractHEOMMatrix, ados::ADOs, Δt::Real, steps::Int; threshold   = 1.0e-6, nonzero_tol = 1.0e-14, verbose::Bool = true, filename::String = "")
```

Using the method based on ODEs : [`OrdinaryDiffEq.jl`](https://github.com/SciML/OrdinaryDiffEq.jl)
```@docs
evolution(M::AbstractHEOMMatrix, ρ0, tlist::AbstractVector; solver = DP5(), reltol::Real = 1.0e-6, abstol::Real = 1.0e-8, maxiters::Real = 1e5, save_everystep::Bool=false, verbose::Bool = true, filename::String = "", SOLVEROptions...)
evolution(M::AbstractHEOMMatrix, ados::ADOs, tlist::AbstractVector; solver = DP5(), reltol::Real = 1.0e-6, abstol::Real = 1.0e-8, maxiters::Real = 1e5, save_everystep::Bool=false, verbose::Bool = true, filename::String = "", SOLVEROptions...)
```

## Steady State
Using the method based on [`LinearSolve.jl`](http://linearsolve.sciml.ai/stable/)
```@docs
SteadyState(M::AbstractHEOMMatrix; solver=UMFPACKFactorization(), verbose::Bool=true, SOLVEROptions...)
```

Using the method based on ODEs : [`OrdinaryDiffEq.jl`](https://github.com/SciML/OrdinaryDiffEq.jl)
```@docs
SteadyState( M::AbstractHEOMMatrix, ρ0; solver = FBDF(autodiff=false), reltol::Real = 1.0e-6, abstol::Real = 1.0e-8, maxiters::Real = 1e5, save_everystep::Bool=false, verbose::Bool = true, SOLVEROptions...)
SteadyState( M::AbstractHEOMMatrix, ados::ADOs; solver = FBDF(autodiff=false), reltol::Real = 1.0e-6, abstol::Real = 1.0e-8, maxiters::Real = 1e5, save_everystep::Bool=false, verbose::Bool = true, SOLVEROptions...)
```

## Spectrum
Power Spectral Density (for bosonic type Bath)
```@docs
PSD
```

Density of States (for fermionic type Bath)
```@docs
DOS
```
# [LinearSolve solvers](@id LS-solvers)

In this page, we list several recommended `alg`orithms provided by [LinearSolve.jl](https://docs.sciml.ai/LinearSolve/stable/) for solving [`steadystate`](@ref) and spectrum in hierarchical equations of motion approach.  

Remember to import `LinearSolve.jl`

```julia
using LinearSolve
```

(click [here](https://docs.sciml.ai/LinearSolve/stable/solvers/solvers/) to see the full `alg`orithm list provided by `LinearSolve.jl`)

### A generic GMRES implementation from Krylov (Default algorithm)

```julia
KrylovJL_GMRES(rtol=1e-12, atol=1e-14)
```

### UMFPACKFactorization
This `alg`orithm performs better when there is more structure to the sparsity pattern (depends on the complexity of your system and baths).

```julia
UMFPACKFactorization()
```

### KLUFactorization
This `alg`orithm performs better when there is less structure to the sparsity pattern (depends on the complexity of your system and baths).

```julia
KLUFactorization()
```

### Pardiso
This `alg`orithm is based on Intel openAPI Math Kernel Library (MKL) Pardiso
!!! note "Note"
    Using this `alg`orithm requires adding the package [Pardiso.jl](https://github.com/JuliaSparse/Pardiso.jl), i.e. `using Pardiso`

```julia
using Pardiso
using LinearSolve
MKLPardisoFactorize()
MKLPardisoIterate()
```

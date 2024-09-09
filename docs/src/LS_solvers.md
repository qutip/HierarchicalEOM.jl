# [LinearSolve solvers](@id LS-solvers)

In this page, we list several recommended solvers provided by [LinearSolve.jl](https://docs.sciml.ai/LinearSolve/stable/) for solving SteadyState and spectrum in hierarchical equations of motion approach.  

Remember to import `LinearSolve.jl`

```julia
using LinearSolve
```

(click [here](https://docs.sciml.ai/LinearSolve/stable/solvers/solvers/) to see the full solver list provided by `LinearSolve.jl`)

### UMFPACKFactorization (Default solver)
This solver performs better when there is more structure to the sparsity pattern (depends on the complexity of your system and baths).

```julia
UMFPACKFactorization()
```

### KLUFactorization
This solver performs better when there is less structure to the sparsity pattern (depends on the complexity of your system and baths).

```julia
KLUFactorization()
```

### A generic BICGSTAB implementation from Krylov

```julia
KrylovJL_BICGSTAB()
```

### Pardiso
This solver is based on Intel openAPI Math Kernel Library (MKL) Pardiso
!!! note "Note"
    Using this solver requires adding the package [Pardiso.jl](https://github.com/JuliaSparse/Pardiso.jl), i.e. `using Pardiso`

```julia
using Pardiso
using LinearSolve
MKLPardisoFactorize()
MKLPardisoIterate()
```

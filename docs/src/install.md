# Installation

## HierarchicalEOM.jl
To install `HierarchicalEOM.jl`, run the following commands inside Julia's interactive session (also known as REPL):
```julia
using Pkg
Pkg.add("HierarchicalEOM")
```
Alternatively, this can also be done in Julia's [Pkg REPL](https://julialang.github.io/Pkg.jl/v1/getting-started/) by pressing the key `]` in the REPL to use the package mode, and then type the following command:
```julia-REPL
(1.8) pkg> add HierarchicalEOM
```
More information about `Julia`'s package manager can be found at [`Pkg.jl`](https://julialang.github.io/Pkg.jl/v1/).  
!!! note "Julia 1.8"
    `HierarchicalEOM.jl` requires Julia 1.8 or higher

To load the package and check the version information, use the command:
```julia
julia> using HierarchicalEOM
julia> HierarchicalEOM.versioninfo()
```

## Other Useful Packages
In order to get a better experience and take full advantage of `HierarchicalEOM`, we recommend to install the following external packages:
### [QuantumOptics.jl](https://qojulia.org/)
`QuantumOptics` framework provides many useful functions to create arbitrary quantum states and operators which can be combined in all the expected ways. Constructing arbitrary Hamiltonians and specifying initial states is also straightforward.  
!!! note "Note" 
    The inputs in `HierarchicalEOM` should be child-type of the standard `Base.AbstractVector` (for vector-type quantum states) and `Base.AbstractMatrix` (for matrix-type quantum operators). Users can still construct the vectors or matrices by standard method and take it as the input of `HierarchicalEOM`.  
    For the users who constructs the quantum objects by `QuantumOptics` package, the standard type of vectors or matrices are stored in the `.data` field of the objects, i.e., `op.data` (where `op` is a `QuantumOptics`-type object). Thus, users should take it as the inputs of `HierarchicalEOM` for objects in `QuantumOptics` package.

!!! compat "Extension for QuantumOptics.jl"
    `HierarchicalEOM.jl` provides an extension to support `QuantumOptics`-type object, but this feature requires `Julia 1.9+` and `HierarchicalEOM 0.3+`. See [here](@ref doc-ext-QuantumOptics) for more details.

### [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/)
`DifferentialEquations` is needed to provide the low-level ODE solvers especially for solving [time evolution](@ref doc-Time-Evolution). For [low dependency usage](https://diffeq.sciml.ai/stable/features/low_dep/), users can use [`OrdinaryDiffEq.jl`](https://github.com/JuliaDiffEq/OrdinaryDiffEq.jl) instead.

### [LinearSolve.jl](http://linearsolve.sciml.ai/stable/)
`LinearSolve` is a unified interface for the linear solving packages of Julia. It interfaces with other packages of the Julia ecosystem to make it easier to test alternative solver packages and pass small types to control algorithm swapping. It is needed to provide the solvers especially for solving [stationary state](@ref doc-Stationary-State) and [spectra](@ref doc-Spectrum) for both bosonic and fermionic systems.

### [JLD2.jl](https://juliaio.github.io/JLD2.jl/stable/)
`JLD2` saves and loads Julia data structures in a format comprising a subset of HDF5. Because the size of matrix in `HierarchicalEOM` is usually super large and leads to long time calculation, we support the functionality for saving and loading the `HierarchicalEOM`-type objects into files by `JLD2 >= 0.4.23`.

### [PyPlot.jl](https://github.com/JuliaPy/PyPlot.jl)
`PyPlot.jl` provides a Julia interface to the [`Matplotlib`](https://matplotlib.org/) plotting library from `Python`, and specifically to the `matplotlib.pyplot` module.
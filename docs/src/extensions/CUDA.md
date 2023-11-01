# [Extension for CUDA.jl](@id doc-ext-CUDA)

This is an extension to support GPU ([`CUDA.jl`](https://github.com/JuliaGPU/CUDA.jl)) acceleration for solving the [time evolution](@ref doc-Time-Evolution) and [spectrum](@ref doc-Spectrum). This improves the execution time and memory usage especially when the HEOMLS matrix is super large.

!!! compat "Compat"
    The described feature requires `Julia 1.9+`.

The functions [`evolution`](@ref doc-Time-Evolution) (only supports ODE method with time-independent system Hamiltonian) and [`spectrum`](@ref doc-Spectrum) will automatically choose to solve on CPU or GPU depend on the type of the sparse matrix in `M::AbstractHEOMLSMatrix` objects (i.e., the type of the field `M.data`). 

```julia
typeof(M.data) <:   SparseMatrixCSC # solve on CPU
typeof(M.data) <: CuSparseMatrixCSC # solve on GPU
```

Therefore, we wrapped several functions in `CUDA` and `CUDA.CUSPARSE` in order to return a new HEOMLS-matrix-type object with `M.data` is in the type of `CuSparseMatrix`, and also change the element type into `ComplexF32` and `Int32` (since GPU performs better in this type). The functions are listed as follows:
- `cu(M::AbstractHEOMLSMatrix)` : Translate `M.data` into the type `CuSparseMatrixCSC{ComplexF32, Int32}`
- `CuSparseMatrixCSC(M::AbstractHEOMLSMatrix)` : Translate `M.data` into the type `CuSparseMatrixCSC{ComplexF32, Int32}`

### Demonstration

The extension will be automatically loaded if user imports the package `CUDA.jl` :

```julia
using CUDA
using HierarchicalEOM
using LinearSolve # to change the solver for better GPU performance
```

### Setup

Here, we demonstrate this extension by using the example of [the single-impurity Anderson model](@ref exp-SIAM). 

```julia
ϵ  = -5
U  = 10
Γ  = 2
μ  = 0
W  = 10
kT = 0.5
N  = 5
tier  = 3

tlist = 0:0.1:10
ωlist = -10:1:10

σm = [0 1; 0  0]
σz = [1 0; 0 -1]
II = [1 0; 0  1]
d_up = kron(     σm, II)
d_dn = kron(-1 * σz, σm)
ρ0   = kron([1 0; 0 0], [1 0; 0 0])
Hsys = ϵ * (d_up' * d_up + d_dn' * d_dn) + U * (d_up' * d_up * d_dn' * d_dn)

bath_up = Fermion_Lorentz_Pade(d_up, Γ, μ, W, kT, N)
bath_dn = Fermion_Lorentz_Pade(d_dn, Γ, μ, W, kT, N)
bath_list = [bath_up, bath_dn]

# even HEOMLS matrix
M_even_cpu = M_Fermion(Hsys, tier, bath_list; verbose=false)
M_even_gpu = cu(M_even_cpu)

# odd HEOMLS matrix
M_odd_cpu  = M_Fermion(Hsys, tier, bath_list, ODD; verbose=false)
M_odd_gpu  = cu(M_odd_cpu)

# solve steady state with CPU
ados_ss = SteadyState(M_even_cpu);
```

!!! note "Note"
    This extension does not support for solving [`SteadyState`](@ref doc-Stationary-State) on GPU since it is not efficient and might get wrong solutions. If you really want to obtain the stationary state with GPU, you can repeatedly solve the [`evolution`](@ref doc-Time-Evolution) until you find it.

### Solving time evolution with CPU

```julia
ados_list_cpu = evolution(M_even_cpu, ρ0, tlist; verbose=false)
```

### Solving time evolution with GPU

```julia
ados_list_gpu = evolution(M_even_gpu, ρ0, tlist; verbose=false)
```

### Solving Spectrum with CPU

```julia
dos_cpu = spectrum(M_odd_cpu, ados_ss, d_up, ωlist; verbose=false)
```

### Solving Spectrum with GPU

```julia
dos_gpu = spectrum(M_odd_gpu, ados_ss, d_up, ωlist; solver=KrylovJL_BICGSTAB(rtol=1f-10, atol=1f-12), verbose=false)
```
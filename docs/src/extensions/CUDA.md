# [Extension for CUDA.jl](@id doc-ext-CUDA)

This is an extension to support GPU ([`CUDA.jl`](https://github.com/JuliaGPU/CUDA.jl)) acceleration for solving the [time evolution](@ref doc-Time-Evolution) and [spectra](@ref doc-Spectrum). This improves the execution time and memory usage especially when the HEOMLS matrix is super large.

The functions of calculating [time evolution](@ref doc-Time-Evolution) (only supports ODE method with time-independent system Hamiltonian) and [spectra](@ref doc-Spectrum) will automatically choose to solve on CPU or GPU depend on the type of the sparse matrix in `M::AbstractHEOMLSMatrix` objects (i.e., the type of the field `M.data`). 

```julia
typeof(M.data) <: SparseMatrixCSC # solve on CPU
typeof(M.data) <: Union{CuSparseMatrixCSC,CuSparseMatrixCSR} # solve on GPU
```

Therefore, we wrapped several functions in `CUDA` and `CUDA.CUSPARSE` in order to not only converting a HEOMLS-matrix-type object into GPU arrays, but also changing the element type and word size (`32` and `64`) since some of the GPUs perform better in `32`-bit. The functions are listed as follows (where input `M` is a `AbstractHEOMLSMatrix`):
- `cu(M, word_size=64)` : Transform `M.data` into `CUDA` sparse arrays with specified `word_size` (default to `64`).
- `CuSparseMatrixCSC(M)` : Translate `M.data` into the type `CuSparseMatrixCSC{eltype(M), Int32}`
- `CuSparseMatrixCSC{T}(M)` : Translate `M.data` into the type `CuSparseMatrixCSC{T, Int32}`
- `CuSparseMatrixCSR(M)` : Translate `M.data` into the type `CuSparseMatrixCSR{eltype(M)}, Int32}`
- `CuSparseMatrixCSR{T}(M)` : Translate `M.data` into the type `CuSparseMatrixCSR{T, Int32}`

### Demonstration

The extension will be automatically loaded if user imports the package `CUDA.jl` :

```julia
using HierarchicalEOM
using CUDA
CUDA.allowscalar(false) # Avoid unexpected scalar indexing
```

### Setup

Here, we demonstrate this extension by using the example of [the single-impurity Anderson model](https://qutip.org/qutip-julia-tutorials/HierarchicalEOM.jl/SIAM.html).

```julia
ϵ  = -5
U  = 10
Γ  = 2
μ  = 0
W  = 10
kT = 0.5
N  = 5
tier = 3

tlist = 0:0.1:10
ωlist = -10:1:10

σm = sigmam()
σz = sigmaz()
II = qeye(2)
d_up = tensor(     σm, II)
d_dn = tensor(-1 * σz, σm)
ψ0   = tensor(basis(2, 0), basis(2, 0))
Hsys = ϵ * (d_up' * d_up + d_dn' * d_dn) + U * (d_up' * d_up * d_dn' * d_dn)

bath_up = Fermion_Lorentz_Pade(d_up, Γ, μ, W, kT, N)
bath_dn = Fermion_Lorentz_Pade(d_dn, Γ, μ, W, kT, N)
bath_list = [bath_up, bath_dn]

# even HEOMLS matrix
M_even_cpu = M_Fermion(Hsys, tier, bath_list)

# odd HEOMLS matrix
M_odd_cpu  = M_Fermion(Hsys, tier, bath_list, ODD)
```

### Generate HEOMLS matrix for GPU computation

```julia
M_even_gpu = cu(M_even_cpu)
M_odd_gpu  = cu(M_odd_cpu)
```

### Solving time evolution with CPU

```julia
ados_list = HEOMsolve(M_even_cpu, ψ0, tlist)
```

### Solving time evolution with GPU

```julia
ados_list = HEOMsolve(M_even_gpu, ψ0, tlist)
```

### Solving steady state with CPU using linear-solve method

```julia
ados_ss = steadystate(M_even_cpu);
```

### Solving steady state with GPU using linear-solve method

```julia
ados_ss = steadystate(M_even_gpu);
```

### Solving steady state with CPU using ODE method

```julia
ados_ss = steadystate(M_even_cpu, ψ0);
```

### Solving steady state with GPU using ODE method

```julia
ados_ss = steadystate(M_even_gpu, ψ0);
```

### Solving spectrum with CPU

```julia
dos = DensityOfStates(M_odd_cpu, ados_ss, d_up, ωlist)
```

### Solving spectrum with GPU

```julia
dos = DensityOfStates(M_odd_gpu, ados_ss, d_up, ωlist)
```

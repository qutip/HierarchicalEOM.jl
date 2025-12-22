# [Lazy Operator Memory Optimization](@id doc-Lazy-Operators)

## Introduction

`HierarchicalEOM.jl` provides a memory-efficient lazy evaluation feature for HEOM Liouvillian superoperator (HEOMLS) matrices using [`TensorProductOperator`](https://docs.sciml.ai/SciMLOperators/stable/) from [`SciMLOperators.jl`](https://github.com/SciML/SciMLOperators.jl). This feature can significantly reduce memory usage for large systems, with savings proportional to the number of auxiliary density operators (ADOs), ``N_{\text{ADO}}``.

## Memory Scaling

### Standard Full Matrix Assembly

By default, HEOMLS matrices are assembled as single large sparse matrices with memory scaling as:
```math
\mathcal{O}(N_{\text{ADO}}^2 \cdot d^4)
```
where ``d`` is the system dimension and ``N_{\text{ADO}}`` grows rapidly with hierarchy tier and number of bath terms.

### Lazy Tensor Product Representation

Instead of assembling a single large matrix, the lazy operator feature represents the HEOM Liouvillian as a sum of Kronecker products:
```math
\mathcal{M} = \sum_i A_i \otimes B_i
```

This reduces memory to:
```math
\mathcal{O}(N_{\text{ADO}} \cdot \text{nnz}(A_i) + d^4 \cdot k)
```
where ``k`` is the number of bath terms. The key insight is that the ``B_i`` matrices (system superoperators) are small (size ``d^2 \times d^2``) and shared across terms, while only the ``A_i`` matrices (ADO connectivity patterns) need to scale with ``N_{\text{ADO}}``.

**Memory savings scale with the number of ADOs**, making this especially beneficial for:
- High hierarchy tiers
- Multiple bath terms
- Large system dimensions
- Long-time dynamics requiring high accuracy

## Mathematical Background: The Lazy Operator Trick

The computational efficiency of lazy operators comes from a mathematical property of Kronecker products. Consider a single term in the HEOM Liouvillian:

```math
\mathcal{M}_i = A_i \otimes B_i
```

where ``A_i`` is an ``N_{\text{ADO}} \times N_{\text{ADO}}`` matrix representing ADO connectivity, and ``B_i`` is a ``d^2 \times d^2`` matrix representing system superoperators.

### Efficient Matrix-Vector Product

When computing the product ``\mathcal{M}_i \vec{v}`` with a vector ``\vec{v}`` of length ``N_{\text{ADO}} \cdot d^2``, we can reshape ``\vec{v}`` into a matrix ``V`` of size ``N_{\text{ADO}} \times d^2`` and use the identity:

```math
(A \otimes B) \mathrm{vec}(V) = \mathrm{vec}(A V B^T)
```

This transforms a single large matrix-vector product into **two smaller matrix-matrix products**:

```math
\mathcal{M}_i \vec{v} = \mathrm{vec}(A_i V B_i^T)
```

### Computational Advantages

1. **Reduced Memory**: Instead of storing the full Kronecker product ``A_i \otimes B_i`` (size ``(N_{\text{ADO}} d^2) \times (N_{\text{ADO}} d^2)``), we only store ``A_i`` (size ``N_{\text{ADO}} \times N_{\text{ADO}}``) and ``B_i`` (size ``d^2 \times d^2``).

2. **Better Parallelization**: Matrix-matrix multiplication (``A_i V`` and then result times ``B_i^T``) can leverage highly optimized BLAS routines (like `gemm`) that are better parallelized than sparse matrix-vector products.

3. **Cache Efficiency**: The two sequential matrix-matrix products have better cache locality than a single sparse matrix-vector product with a very large matrix.

### Full HEOM Liouvillian

For the complete HEOM Liouvillian with multiple terms:

```math
\mathcal{M} \vec{v} = \sum_i (A_i \otimes B_i) \vec{v} = \sum_i \mathrm{vec}(A_i V B_i^T)
```

The `Val(:combine)` mode groups terms with identical ``B_i`` to further reduce the number of operations and memory overhead.

## The `assemble` Parameter

The `assemble` keyword parameter controls how the HEOMLS matrix is constructed. It is available in:
- [`M_Boson`](@ref)
- [`M_Fermion`](@ref)
- [`M_Boson_Fermion`](@ref)

### Three Assembly Modes

#### `Val(:full)` - Full Sparse Matrix (Default)
Assembles the complete HEOMLS matrix as a single sparse matrix. This is the default behavior for backward compatibility.

```julia
L_full = M_Boson(Hsys, tier, bath; assemble = Val(:full))
```

#### `Val(:combine)` - Combined Lazy Operators (Recommended)
Combines terms with identical system operators (``B_i``) but does not materialize the full matrix. This provides excellent memory savings while maintaining computational efficiency. The lazy representation can leverage matrix-matrix multiplication patterns, which are more parallelizable than the sparse matrix-vector products used by full assembly.

```julia
L_lazy = M_Boson(Hsys, tier, bath; assemble = Val(:combine))
```

**Recommended for most applications, especially memory-constrained systems.**

#### `Val(:none)` - Individual Lazy Operators (For Sparsity Analysis)
Keeps all tensor product terms separate without any combination. This mode actually uses **more memory** than `Val(:combine)` because `Val(:combine)` reduces the number of terms in the `AddedOperator` by grouping terms with identical system operators. However, `Val(:none)` provides flexibility for investigating sparsity patterns and analyzing the structure of individual tensor product terms.

```julia
L_lazy_uncombined = M_Boson(Hsys, tier, bath; assemble = Val(:none))
```

!!! note "Memory Usage: `:none` vs `:combine`"
    Contrary to what the name might suggest, `Val(:none)` does not provide the maximum memory savings. In fact, `Val(:combine)` is more memory-efficient because it reduces the number of operator terms that need to be stored. Use `Val(:none)` primarily for flexibility in analyzing the operator structure, not for memory optimization.

## Combining with Importance Threshold

The lazy operator feature is **completely compatible** with the [importance threshold](@ref doc-Importance-Value-and-Threshold) option, allowing you to combine two powerful memory optimization techniques:

1. **Lazy operators** reduce memory by avoiding full matrix assembly (saves proportional to ``N_{\text{ADO}}``)
2. **Importance threshold** reduces the number of ADOs by filtering out less important ones (reduces ``N_{\text{ADO}}`` itself)

These optimizations are complementary and can be used together:

```julia
using HierarchicalEOM

# System setup
Hsys = 0.25 * sigmaz() + 0.5 * sigmax()
bath = Boson_DrudeLorentz_Pade(sigmaz(), 0.01, 0.5, 0.5, 3)

# Combine lazy operators with importance filtering
# threshold = 1e-6 filters out ADOs with importance < 1e-6
L_optimized = M_Boson(Hsys, 6, bath; 
                      assemble = Val(:combine),  # Lazy operators
                      threshold = 1e-6)          # Importance filtering

println("Number of ADOs: $(L_optimized.N)")
println("Memory usage: $(round(Base.summarysize(L_optimized.data) / 1024^2, digits=2)) MB")

# Compare with full assembly and no filtering
L_full = M_Boson(Hsys, 6, bath; assemble = Val(:full), threshold = 0.0)
println("Full matrix ADOs: $(L_full.N)")
println("Full matrix memory: $(round(Base.summarysize(L_full.data) / 1024^2, digits=2)) MB")
```

The actual memory savings depend on the specific system, hierarchy tier, and threshold value. Both techniques work together to reduce overall memory consumption.

!!! tip "Recommended for Large Systems"
    For systems with high hierarchy tiers (≥6) or multiple baths, we recommend:
    ```julia
    L = M_Boson(Hsys, tier, bath; 
                assemble = Val(:combine),
                threshold = 1e-6)  # Adjust threshold based on required accuracy
    ```
    This combination provides excellent memory efficiency while maintaining accuracy.

## Usage Examples

### Basic Setup

```julia
using HierarchicalEOM

# System Hamiltonian (two-level system)
Hsys = 0.25 * sigmaz() + 0.5 * sigmax()

# Bosonic bath with Drude-Lorentz spectral density
bath = Boson_DrudeLorentz_Pade(sigmaz(), 0.01, 0.5, 0.5, 3)

# Standard full matrix (default)
L_full = M_Boson(Hsys, 5, bath)

# Memory-efficient lazy matrix (recommended)
L_lazy = M_Boson(Hsys, 5, bath; assemble = Val(:combine))

# Uncombined lazy operators (for sparsity analysis)
L_lazy_uncombined = M_Boson(Hsys, 5, bath; assemble = Val(:none))
```

### Demonstrating Memory Savings

Use `Base.summarysize` to measure and compare memory usage:

```julia
using HierarchicalEOM

# Example system
Hsys = 0.25 * sigmaz() + 0.5 * sigmax()
bath = Boson_DrudeLorentz_Pade(sigmaz(), 0.01, 0.5, 0.5, 3)

# Create matrices with different assembly modes
L_full = M_Boson(Hsys, 6, bath; assemble = Val(:full))
L_combine = M_Boson(Hsys, 6, bath; assemble = Val(:combine))
L_none = M_Boson(Hsys, 6, bath; assemble = Val(:none))

# Measure memory usage
mem_full = Base.summarysize(L_full.data) / 1024^2  # Convert to MB
mem_combine = Base.summarysize(L_combine.data) / 1024^2
mem_none = Base.summarysize(L_none.data) / 1024^2

println("Memory usage:")
println("  Full assembly:      $(round(mem_full, digits=2)) MB")
println("  Combined lazy:      $(round(mem_combine, digits=2)) MB")
println("  Individual lazy:    $(round(mem_none, digits=2)) MB")
println("Memory reduction:")
println("  Combined vs Full:   $(round(100 * (1 - mem_combine/mem_full), digits=1))%")
println("  Individual vs Full: $(round(100 * (1 - mem_none/mem_full), digits=1))%")
```

For a tier-6 system with 84 ADOs, you might see memory reductions of 60-80% with lazy operators. **The savings increase dramatically with higher tiers and more ADOs.**

### Larger System Example

```julia
# System with higher tier and multiple baths
Hsys = 0.25 * sigmaz() + 0.5 * sigmax()
bath1 = Boson_DrudeLorentz_Pade(sigmaz(), 0.01, 0.5, 0.5, 3)
bath2 = Boson_DrudeLorentz_Pade(sigmax(), 0.01, 0.3, 0.5, 3)

# High tier for better accuracy
tier = 8

# Full assembly (will use significant memory)
L_full = M_Boson(Hsys, tier, [bath1, bath2]; assemble = Val(:full))

# Lazy assembly (recommended for large systems)
L_lazy = M_Boson(Hsys, tier, [bath1, bath2]; assemble = Val(:combine))

println("ADOs: $(L_full.N)")
println("Full matrix memory: $(round(Base.summarysize(L_full.data) / 1024^2, digits=2)) MB")
println("Lazy matrix memory: $(round(Base.summarysize(L_lazy.data) / 1024^2, digits=2)) MB")
```

## Time Evolution with Lazy Operators

Lazy operators work seamlessly with time evolution in `HierarchicalEOM.jl`. No code changes are needed:

```julia
using HierarchicalEOM

# Setup
Hsys = 0.25 * sigmaz() + 0.5 * sigmax()
bath = Boson_DrudeLorentz_Pade(sigmaz(), 0.01, 0.5, 0.5, 3)
L_lazy = M_Boson(Hsys, 6, bath; assemble = Val(:combine))

# Initial state
ψ0 = basis(2, 0)  # Ground state
ρ0 = ψ0 * ψ0'

# Time evolution (works directly with lazy operators)
tlist = 0:0.1:10
sol = HEOMsolve(L_lazy, ρ0, tlist; e_ops = [sigmaz(), sigmax()])

# Access results as usual
expect_z = sol.expect[1, :]
expect_x = sol.expect[2, :]
```

See [Time Evolution](@ref doc-Time-Evolution) for more details.

## Converting Lazy to Full Matrix

Sometimes you may need the full materialized matrix, for example when computing eigenvalues or for analysis. Use `SciMLOperators.concretize()` to convert a lazy operator back to a full sparse matrix:

```julia
using HierarchicalEOM
using SciMLOperators: concretize

# Create lazy HEOMLS matrix
Hsys = 0.25 * sigmaz() + 0.5 * sigmax()
bath = Boson_DrudeLorentz_Pade(sigmaz(), 0.01, 0.5, 0.5, 3)
L_lazy = M_Boson(Hsys, 5, bath; assemble = Val(:combine))

# Convert to full sparse matrix when needed
L_materialized = concretize(L_lazy.data)

# Now you can use standard linear algebra operations
# For example, compute eigenvalues
using LinearAlgebra
eigenvalues = eigvals(Matrix(L_materialized))
```

!!! note "When to Use Full vs Lazy"
    - **Use lazy operators** (`Val(:combine)`) when:
      - Memory is constrained
      - System has high tier or many ADOs
      - Performing time evolution or steady-state calculations
      - Working with systems where memory is the bottleneck
    
    - **Use full assembly** (`Val(:full)`) when:
      - Backward compatibility is required
      - Require eigenvalue/eigenvector analysis
      - Working with specialized linear algebra routines that require concrete matrices

## GPU Compatibility

Lazy operators are fully compatible with GPU acceleration via [`CUDA.jl`](https://github.com/JuliaGPU/CUDA.jl). Simply convert the matrix to GPU arrays using `cu()`:

```julia
using HierarchicalEOM
using CUDA

# Create lazy HEOMLS matrix
Hsys = 0.25 * sigmaz() + 0.5 * sigmax()
bath = Boson_DrudeLorentz_Pade(sigmaz(), 0.01, 0.5, 0.5, 3)
L_lazy = M_Boson(Hsys, 5, bath; assemble = Val(:combine))

# Convert to GPU
L_gpu = cu(L_lazy)

# Time evolution on GPU
ψ0 = basis(2, 0)
ρ0 = ψ0 * ψ0'
tlist = 0:0.1:10
sol = HEOMsolve(L_gpu, ρ0, tlist; e_ops = [sigmaz()])
```

See [Extension for CUDA.jl](@ref doc-ext-CUDA) for more details on GPU acceleration.

## Performance Considerations

Lazy operators provide substantial memory savings with competitive computational performance:

### Memory vs Computation

- **Memory**: Lazy operators reduce memory by a factor proportional to ``N_{\text{ADO}}``, often achieving 60-90% reduction for high-tier systems.
- **Computation**: The lazy tensor product representation enables matrix-matrix multiplication patterns, which can be more parallelizable than sparse matrix-vector products. Performance depends on the specific system size and available computational resources.

### Recommendations

1. **For memory-constrained systems**: Always use `assemble = Val(:combine)`. The memory savings are substantial and computational performance remains competitive.

2. **For small systems**: Either `assemble = Val(:full)` or `assemble = Val(:combine)` can be used, as memory is typically not a concern for small systems.

3. **For sparsity analysis**: Use `assemble = Val(:none)` to keep all tensor product terms separate for investigating the operator structure. Note that this uses more memory than `Val(:combine)` due to the larger number of terms stored.

4. **Profile your specific case**: Use `@time` and `Base.summarysize` to measure the trade-offs for your particular system:

```julia
using HierarchicalEOM

Hsys = 0.25 * sigmaz() + 0.5 * sigmax()
bath = Boson_DrudeLorentz_Pade(sigmaz(), 0.01, 0.5, 0.5, 3)
ψ0 = basis(2, 0)
ρ0 = ψ0 * ψ0'
tlist = 0:0.1:5

# Time full assembly
L_full = M_Boson(Hsys, 6, bath; assemble = Val(:full))
@time sol_full = HEOMsolve(L_full, ρ0, tlist; e_ops = [sigmaz()])

# Time lazy assembly
L_lazy = M_Boson(Hsys, 6, bath; assemble = Val(:combine))
@time sol_lazy = HEOMsolve(L_lazy, ρ0, tlist; e_ops = [sigmaz()])

# Compare memory
println("Memory - Full: $(Base.summarysize(L_full.data) / 1024^2) MB")
println("Memory - Lazy: $(Base.summarysize(L_lazy.data) / 1024^2) MB")
```

## Summary

The lazy operator feature provides substantial memory savings for HEOM calculations by representing the Liouvillian as a sum of Kronecker products rather than a single large sparse matrix. The `assemble` parameter gives you control over the representation:

- `Val(:full)`: Default, assembles full sparse matrix (for backward compatibility)
- `Val(:combine)`: **Recommended** for most applications, provides best memory savings with competitive performance
- `Val(:none)`: For sparsity analysis (uses more memory than `:combine` but provides flexibility)

The memory savings scale with ``N_{\text{ADO}}``, making this feature essential for high-tier calculations, multiple baths, or long-time dynamics.

# [Memory Optimization with Lazy Operators](@id doc-Lazy-Operators)

## Introduction

`HierarchicalEOM.jl` provides a memory-efficient lazy evaluation feature for HEOM Liouvillian superoperator (HEOMLS) matrices using [`SciMLOperators.TensorProductOperator`](https://docs.sciml.ai/SciMLOperators/stable/). This feature can significantly reduce memory usage for large systems, with savings proportional to the number of auxiliary density operators (ADOs), ``N_{\text{ADO}}``.

By default, HEOMLS matrices are assembled as single large sparse matrices with memory scaling as:

```math
\mathcal{O}(N_{\text{ADO}}^2 \cdot d^4)
```

where ``d`` is the system dimension and ``N_{\text{ADO}}`` grows rapidly with hierarchy tier and number of bath terms.

Instead of assembling a single large matrix, the lazy operator feature represents the HEOM Liouvillian as a sum of Kronecker products:

```math
\mathcal{M} = \sum_i A_i \otimes B_i
```

This reduces memory to:

```math
\mathcal{O}(N_{\text{ADO}} \cdot \text{nnz}(A_i) + d^4 \cdot k)
```

where ``k`` is the number of bath terms. The key insight is that the ``B_i`` matrices (system superoperators) are small (size ``d^2 \times d^2``) and shared across terms, while only the ``A_i`` matrices (ADO connectivity patterns) need to scale with ``N_{\text{ADO}}``.

Memory savings scale with the number of ADOs, making this especially beneficial for:

- High hierarchy tiers
- Multiple bath terms
- Large system dimensions
- Long-time dynamics requiring high accuracy

## The `assemble` keyword argument

The `assemble` keyword argument controls how the HEOMLS matrix is constructed. It is available in:

- [`M_Boson`](@ref)
- [`M_Fermion`](@ref)
- [`M_Boson_Fermion`](@ref)

#### `assemble = Val(:full)` - Full Sparse Matrix (Default)

Assembles the complete HEOMLS matrix as a single sparse matrix. This is the default behavior for backward compatibility.

```julia
L_full = M_Boson(Hsys, tier, bath; assemble = Val(:full))
```

#### `assemble = Val(:combine)` - Combined Lazy Operators (Recommended)

Combines terms with identical system operators (``B_i``) but does not materialize the full matrix. This provides excellent memory savings while maintaining computational efficiency. The lazy representation can leverage matrix-matrix multiplication patterns, which are more parallelizable than the sparse matrix-vector products used by full assembly.

```julia
L_combine = M_Boson(Hsys, tier, bath; assemble = Val(:combine))
```

Recommended for most applications, especially memory-constrained devices.

#### `assemble = Val(:none)` - Individual Lazy Operators (For Sparsity Analysis)

Keeps all tensor product terms separate without any combination. This mode actually uses more memory than `Val(:combine)` because `Val(:combine)` reduces the number of terms in the `SciMLOperators.AddedOperator` by grouping terms with identical system operators. However, `Val(:none)` provides flexibility for investigating sparsity patterns and analyzing the structure of individual tensor product terms.

```julia
L_none = M_Boson(Hsys, tier, bath; assemble = Val(:none))
```

## Demonstrating Memory Savings

```@setup lazy-operators
using HierarchicalEOM
```

```@example lazy-operators 
# spin-boson model
Hsys = 1.0 * sigmaz() + 0.25 * sigmax()

Δ  = 0.1 # coupling strength
W  = 0.2 # band-width  of the environment
kT = 0.5 # the product of the Boltzmann constant k and the absolute temperature T
N  = 10  # Number of exponential terms
baths = [
    Boson_DrudeLorentz_Pade(sigmax(), 0.01, 0.5, 0.5, 3),
    Boson_DrudeLorentz_Pade(sigmax(), 0.01, 0.5, 0.5, 3),
]

# Create matrices with different assembly modes
tier = 10
L_full    = M_Boson(Hsys, tier, baths; assemble=Val(:full),    verbose=false)
L_combine = M_Boson(Hsys, tier, baths; assemble=Val(:combine), verbose=false)
L_none    = M_Boson(Hsys, tier, baths; assemble=Val(:none),    verbose=false)
nothing # hide
```

Use `Base.summarysize` to measure and compare memory usage (convert to `MB`)

```@example lazy-operators
mem_full = Base.summarysize(L_full.data) / 1024^2
mem_combine = Base.summarysize(L_combine.data) / 1024^2
mem_none = Base.summarysize(L_none.data) / 1024^2

println("Number of ADOs: $(L_full.N)")
println("Memory usage:")
println("  full    : $(round(mem_full, digits=3)) MB")
println("  combine : $(round(mem_combine, digits=3)) MB")
println("  none    : $(round(mem_none, digits=3)) MB")
println("Memory reduction:")
println("  combine vs full : $(round(100 * (1 - mem_combine/mem_full), digits=1))%")
println("   none   vs full : $(round(100 * (1 - mem_none/mem_full), digits=1))%")
```

In this case, you might see memory reductions of 80-90% with lazy operators. The savings increase dramatically with higher tiers and more ADOs.

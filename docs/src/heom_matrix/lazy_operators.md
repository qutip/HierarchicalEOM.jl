# [Memory Optimization with Lazy Operators](@id doc-Lazy-Operators)

## Introduction

`HierarchicalEOM.jl` provides a memory-efficient lazy evaluation feature for HEOM Liouvillian superoperator (HEOMLS) matrices using [`SciMLOperators.TensorProductOperator`](https://docs.sciml.ai/SciMLOperators/stable/). This feature can significantly reduce memory usage for large systems.

By default, a HEOMLS matrix is assembled as a single large sparse matrix with space complexity (memory scaling) as:

```math
\mathcal{O}(N_{\text{ADO}}^2 \cdot d^4)
```

where ``d`` is the system dimension, and ``N_{\text{ADO}}`` is the number of auxiliary density operators (ADOs) which grows rapidly with hierarchy tier and number of bath terms.

Instead of assembling HEOMLS into a single large sparse matrix, the lazy operator feature represents it as a sum of Kronecker products:

```math
\hat{\mathcal{M}} = \sum_i A_i \otimes B_i,
```

where ``B_i`` represent the system coupling superoperators, and ``A_i`` encodes the hierarchy (ADO) connectivity pattern and corresponding prefix values of ``B_i``. In this case, the space complexity (memory allocation) has been reduced to

```math
\mathcal{O}(N_{\text{ADO}}^2 + d^4),
```

which significantly lowers the memory requirements for high-dimensional simulations. The key insight is that the ``B_i`` (system coupling superoperator) matrices are small (size ``d^2 \times d^2``) and shared across different connections (coupling) between ADOs, while only the ``A_i`` matrices (hierarchy connectivity patterns) need to scale with ``N_{\text{ADO}}^2``.

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

### `assemble = Val(:full)` - Full Sparse Matrix (Default)

Assembles the complete HEOMLS matrix as a single sparse matrix. This is the default behavior for backward compatibility.

```julia
L_full = M_Boson(Hsys, tier, bath; assemble = Val(:full))
```

### `assemble = Val(:combine)` - Combined Lazy Operators (Recommended)

Combines terms with identical system coupling superoperators (``B_i``) but does not materialize the full matrix. This provides excellent memory savings while maintaining computational efficiency. The lazy representation can leverage matrix-matrix multiplication patterns, which are more parallelizable than the sparse matrix-vector products used by full assembly.

```julia
L_combine = M_Boson(Hsys, tier, bath; assemble = Val(:combine))
```

Recommended for most applications, especially memory-constrained devices.

### `assemble = Val(:none)` - Individual Lazy Operators (For Sparsity Analysis)

Keeps all tensor product terms separate without any combination. This mode actually uses more memory than `Val(:combine)` because `Val(:combine)` reduces the number of terms in the `SciMLOperators.AddedOperator` by grouping terms with identical system coupling superoperators. However, `Val(:none)` provides flexibility for investigating sparsity patterns and analyzing the structure of individual tensor product terms.

```julia
L_none = M_Boson(Hsys, tier, bath; assemble = Val(:none))
```

## Demonstrating Memory Savings

```@setup lazy-operators
using HierarchicalEOM
using SparseArrays
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

print(
    "System dimension (d)     : $(size(Hsys, 1))\n",
    "Number of ADOs   (N_ADO) : $(L_full.N)\n",
    "L_full matrix size       : $(size(L_full, 1)) × $(size(L_full, 2))\n",
    "L_full non-zero elements : $(nnz(L_full.data.A))\n",
    "\n",
    "Memory usage:\n",
    "  full    : $(round(mem_full, digits=3)) MB\n",
    "  combine : $(round(mem_combine, digits=3)) MB\n",
    "  none    : $(round(mem_none, digits=3)) MB\n",
    "Memory reduction:\n",
    "  combine vs full : $(round(100 * (1 - mem_combine/mem_full), digits=1))%\n",
    "   none   vs full : $(round(100 * (1 - mem_none/mem_full), digits=1))%\n",
)
```

In this case, you might see memory reductions of 80-90% with lazy operators. The savings increase dramatically with higher tiers and more ADOs.

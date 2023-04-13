# [Hierarchy Dictionary](@id doc-Hierarchy-Dictionary)
## Introduction
For hierarchical equations of motions, there are many indices the users have to deal with including the indices of the `Exponent` in [bosonic baths](@ref doc-Bosonic-Bath), the `Exponent` in [fermionic baths](@ref doc-Fermionic-Bath), and the ADOs formed by the hierarchy.

With the auxiliary density operators ``\rho_{\textbf{j}\vert\textbf{q}}^{(m,n,p)}``, we use the following keywords :
 - `idx` : the index of the [auxiliary density operators](@ref doc-ADOs)
 - `lvl` : the level ``m`` (``n``) of the bosonic (fermionic) hierarchy
 - `nvec` : object [`Nvec`](@ref) which stores the number of existence for each multi-index ensemble ``j`` (``q``) in vector ``\textbf{j}`` (``\textbf{q}``).

## Dictionary for Pure Bosonic or Fermionic Baths
An object which contains all dictionaries for pure (bosonic or fermionic) bath-ADOs hierarchy, defined as:

[`struct HierarchyDict <: AbstractHierarchyDict`](@ref HierarchyDict)

[`HierarchyDict`](@ref) can be obtained from the field `.hierarchy` in [`M_Boson`](@ref) or [`M_Fermion`](@ref), and it contains the following fields :
 - `idx2nvec` : Return the `Nvec` from a given index of ADO
 - `nvec2idx` : Return the index of ADO from a given `Nvec`
 - `lvl2idx` : Return the list of ADO-indices from a given hierarchy level
 - `bathPtr` : Records the tuple ``(\alpha, k)`` for each position in `Nvec`, where ``\alpha`` and ``k`` represents the ``k``-th exponential-expansion term of the ``\alpha``-th bath.

```julia
# HEOMLS for bosonic baths
M::M_Boson
HDict = M.hierarchy

# HEOMLS for fermionic baths
M::M_Fermion
HDict = M.hierarchy

# obtain the nvec corresponds to 10-th ADO
nvec = HDict.idx2nvec[10]

# obtain the index of the ADO corresponds to the given nvec
nvec::Nvec
idx = HDict.nvec2idx[nvec]

# obtain a list of indices which corresponds to all ADOs in 3rd-level of hierarchy
idx_list = HDict.lvl2idx[3] 
```

`Heom.jl` also provides a function [`getIndexEnsemble(nvec, bathPtr)`](@ref) to obtain the index of the [`Exponent`](@ref) and it's corresponding index of bath:
```julia
# HEOMLS
M::M_Boson
M::Fermion

HDict = M.hierarchy

# auxiliary density operators
ados::ADOs

for (idx, ado) in enumerate(ados)
    ado # the corresponding auxiliary density operator for idx

    # obtain the nvec corresponds to ado
    nvec = HDict.idx2nvec[idx]

    for (α, k, n) in getIndexEnsemble(nvec, HDict.bathPtr)
        α  # index of the bath
        k  # the index of the exponential-expansion term in α-th bath
        n  # the repetition number of the ensemble {α, k} in vector j (or q) in ADOs
        exponent = M.bath[α][k]  # the k-th exponential-expansion term in α-th bath

        # do some calculations you want
    end
end
```

## Dictionary for Mixed Bosonic and Fermionic Baths
An object which contains all dictionaries for mixed (bosonic and fermionic) bath-ADOs hierarchy, defined as:

[`struct MixHierarchyDict <: AbstractHierarchyDict`](@ref MixHierarchyDict)

[`MixHierarchyDict`](@ref) can be obtained from the field `.hierarchy` in [`M_Boson_Fermion`](@ref), and it contains the following fields :
 - `idx2nvec` : Return the tuple `(Nvec_b, Nvec_f)` from a given index of ADO, where `b` represents boson and `f` represents fermion
 - `nvec2idx` : Return the index from a given tuple `(Nvec_b, Nvec_f)`, where `b` represents boson and `f` represents fermion
 - `Blvl2idx` : Return the list of ADO-indices from a given bosonic-hierarchy level
 - `Flvl2idx` : Return the list of ADO-indices from a given fermionic-hierarchy level
 - `bosonPtr` : Records the tuple ``(\alpha, k)`` for each position in `Nvec_b`, where ``\alpha`` and ``k`` represents the ``k``-th exponential-expansion term of the ``\alpha``-th bosonic bath.
 - `fermionPtr` : Records the tuple ``(\alpha, k)`` for each position in `Nvec_f`, where ``\alpha`` and ``k`` represents the ``k``-th exponential-expansion term of the ``\alpha``-th fermionic bath.

```julia
# HEOMLS 
M::M_Boson_Fermion
HDict = M.hierarchy

# obtain the nvec(s) correspond to 10-th ADO
nvec_b, nvec_f = HDict.idx2nvec[10]

# obtain the index of the ADO corresponds to the given nvec
nvec_b::Nvec
nvec_f::Nvec
idx = HDict.nvec2idx[(nvec_b, nvec_f)]

# obtain a list of indices which corresponds to all ADOs in 3rd-bosonic-level of hierarchy
idx_list = HDict.Blvl2idx[3] 

# obtain a list of indices which corresponds to all ADOs in 4rd-fermionic-level of hierarchy
idx_list = HDict.Flvl2idx[4] 
```

`Heom.jl` also provides a function [`getIndexEnsemble(nvec, bathPtr)`](@ref) to obtain the index of the [`Exponent`](@ref) and it's corresponding index of bath:
```julia
# HEOMLS
M::M_Boson_Fermion

HDict = M.hierarchy

# auxiliary density operators
ados::ADOs

for (idx, ado) in enumerate(ados)
    ado # the corresponding auxiliary density operator for idx

    # obtain the nvec(s) correspond to ado
    nvec_b, nvec_f = HDict.idx2nvec[idx]

    # bosonic bath indices
    for (β, k, n) in getIndexEnsemble(nvec_b, HDict.bosonPtr)
        β  # index of the bosonic bath
        k  # the index of the exponential-expansion term in β-th bosonic bath
        nb # the repetition number of the ensemble {β, k} in vector j in ADOs
        exponent = M.Bbath[β][k]  # the k-th exponential-expansion term in β-th bosonic bath

        # do some calculations you want
    end

    # fermionic bath indices
    for (α, h, n) in getIndexEnsemble(nvec_f, HDict.fermionPtr)
        α  # index of the fermionic bath
        h  # the index of the exponential-expansion term in α-th fermionic bath
        nf # the repetition number of the ensemble {α, h} in vector q in ADOs
        exponent = M.Bbath[α][h]  # the h-th exponential-expansion term in α-th fermionic bath

        # do some calculations you want
    end
end
```
"""
    struct HierarchyDict
An object which contains all dictionaries for bath-ADOs hierarchy.

The `ado` (or `n_vector`) denotes a set of integers:
```math
\\{ n_{11}, ..., n_{\\nu k}, ... \\}
```
where ``n_{\\nu k} \\geq 0`` associated with the ``k``-th exponential-expansion term in the ``\\nu``-th bath.

The hierarchy level (``L``) for an `ado` is given by ``L=\\sum_{\\nu, k} n_{\\nu k}``

# Fields
- `idx2ado` : Return the `ado` (`n_vector`) from a given index
- `ado2idx` : Return the index from a given `ado` (`n_vector`)
- `lvl2idx` : Return the list of indices from a given level
- `bathPtr` : The indices for the exponential-expansion terms of the `N`-th bath in `ado` (`n_vector`) can be obtained by `bathPtr[N]:(bathPtr[N+1] - 1)`
"""
struct HierarchyDict
    idx2ado::Vector{Vector{Int}}
    ado2idx::Dict{Vector{Int}, Int}
    lvl2idx::Dict{Int, Vector{Int}}
    bathPtr::Vector{Int}
end

# generate index to ado vector
function _IDX2ADO(n_vec::Vector{Int}, N_exc::Int)
    len = length(n_vec)
    ado = zeros(Int, len)
    result = [copy(ado)]
    nexc = 0

    while true
        idx = len
        ado[end] += 1
        nexc += 1
        if ado[idx] < n_vec[idx]
            push!(result, copy(ado))
        end
        while (nexc == N_exc) || (ado[idx] == n_vec[idx])
            #ado[idx] = 0
            idx -= 1
            if idx < 1
                return result
            end

            nexc -= ado[idx + 1] - 1
            ado[idx + 1] = 0
            ado[idx] += 1
            if ado[idx] < n_vec[idx]
                push!(result, copy(ado))
            end
        end
    end
end

function _genHierarchyDict(n_vec::Vector{Int}, N_exc::Int)
    idx2ado = _IDX2ADO(n_vec, N_exc)
    ado2idx = Dict{Vector{Int}, Int}()
    
    # create lvl2idx
    lvl2idx = Dict{Int, Vector{Int}}()
    for level in 0:N_exc
        lvl2idx[level] = []
    end

    for (idx, ado) in enumerate(idx2ado)
        level = sum(ado)
        push!(lvl2idx[level], idx)
        ado2idx[ado] = idx        
    end

    return length(idx2ado), idx2ado, ado2idx, lvl2idx
end

function genBathHierarchy(B::Vector{T}, tier::Int, dim::Int) where T <: AbstractBath
    Nterm   = 0
    bathPtr = [1]

    if T == BosonBath
        baths = AbstractBosonBath[]
        for b in B
            if b.dim != dim 
                error("The matrix size of the bosonic bath coupling operators are not consistent.")
            end
            push!(baths, b.bath...)
            Nterm += b.Nterm
            push!(bathPtr, Nterm + 1)
        end
        n_vec = fill((tier + 1), Nterm)
    
    elseif T == FermionBath
        baths = AbstractFermionBath[]
        for b in B
            if b.dim != dim 
                error("The matrix size of the fermionic bath coupling operators are not consistent.")
            end
            push!(baths, b.bath...)
            Nterm += b.Nterm
            push!(bathPtr, Nterm + 1)
        end
        n_vec = fill(2, Nterm)
    end

    Nado, idx2ado, ado2idx, lvl2idx = _genHierarchyDict(n_vec, tier)
    hierarchy = HierarchyDict(idx2ado, ado2idx, lvl2idx, bathPtr)

    return Nado, baths, hierarchy
end
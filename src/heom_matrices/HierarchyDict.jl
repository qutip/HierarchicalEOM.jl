"""
    struct HierarchyDict
An object which contains all dictionaries for bath-ADOs hierarchy.

The `nvec` (``\\vec{n}``) denotes a set of integers:
```math
\\{ n_{11}, ..., n_{\\nu k}, ... \\}
```
where ``n_{\\nu k} \\geq 0`` associated with the ``k``-th exponential-expansion term in the ``\\nu``-th bath.

The hierarchy level (``L``) for an `nvec` is given by ``L=\\sum_{\\nu, k} n_{\\nu k}``

# Fields
- `idx2nvec` : Return the `nvec` from a given index
- `nvec2idx` : Return the index from a given `nvec`
- `lvl2idx` : Return the list of indices from a given level
- `bathPtr` : Records the tuple ``(k, \\nu)`` for each position in `n_vector`, where ``k`` and ``\\nu`` represents the ``\\nu``-th exponential-expansion term of the ``k``-th bath.
"""
struct HierarchyDict
    idx2nvec::Vector{Vector{Int}}
    nvec2idx::Dict{Vector{Int}, Int}
    lvl2idx::Dict{Int, Vector{Int}}
    bathPtr::Vector{Tuple}
end

# generate index to n vector
function _Idx2Nvec(n_vec::Vector{Int}, N_exc::Int)
    len = length(n_vec)
    nvec = zeros(Int, len)
    result = [copy(nvec)]
    nexc = 0

    while true
        idx = len
        nvec[end] += 1
        nexc += 1
        if nvec[idx] < n_vec[idx]
            push!(result, copy(nvec))
        end
        while (nexc == N_exc) || (nvec[idx] == n_vec[idx])
            #nvec[idx] = 0
            idx -= 1
            if idx < 1
                return result
            end

            nexc -= nvec[idx + 1] - 1
            nvec[idx + 1] = 0
            nvec[idx] += 1
            if nvec[idx] < n_vec[idx]
                push!(result, copy(nvec))
            end
        end
    end
end

function genBathHierarchy(B::Vector{T}, tier::Int, dim::Int; threshold::Real=0.0) where T <: AbstractBath
    Nterm   = 0
    bathPtr = Tuple[]

    if T == BosonBath
        baths = AbstractBosonBath[]
        for (k, b) in enumerate(B)
            if b.dim != dim 
                error("The matrix size of the bosonic bath coupling operators are not consistent.")
            end
            push!(baths, b.bath...)
            for ν in 1:b.Nterm
                push!(bathPtr, (k, ν))
            end
            Nterm += b.Nterm
        end
        n_vec = fill((tier + 1), Nterm)
    
    elseif T == FermionBath
        baths = AbstractFermionBath[]
        for (k, b) in enumerate(B)
            if b.dim != dim 
                error("The matrix size of the fermionic bath coupling operators are not consistent.")
            end
            push!(baths, b.bath...)
            for ν in 1:b.Nterm
                push!(bathPtr, (k, ν))
            end
            Nterm += b.Nterm
        end
        if tier == 0
            n_vec = fill(1, Nterm)
        elseif tier >= 1
            n_vec = fill(2, Nterm)
        end
    end

    # create idx2nvec
    idx2nvec = _Idx2Nvec(n_vec, tier)

    # create lvl2idx and nvec2idx
    lvl2idx = Dict{Int, Vector{Int}}()
    nvec2idx = Dict{Vector{Int}, Int}()
    for level in 0:tier
        lvl2idx[level] = []
    end
    for (idx, nvec) in enumerate(idx2nvec)
        level = sum(nvec)
        push!(lvl2idx[level], idx)
        nvec2idx[nvec] = idx
    end

    hierarchy = HierarchyDict(idx2nvec, nvec2idx, lvl2idx, bathPtr)

    return length(idx2nvec), baths, hierarchy
end
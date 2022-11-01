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
function _Idx2Nvec(n_max::Vector{Int}, N_exc::Int)
    len = length(n_max)
    nvec = zeros(Int, len)
    result = [copy(nvec)]
    nexc = 0

    while true
        idx = len
        nvec[end] += 1
        nexc += 1
        if nvec[idx] < n_max[idx]
            push!(result, copy(nvec))
        end
        while (nexc == N_exc) || (nvec[idx] == n_max[idx])
            #nvec[idx] = 0
            idx -= 1
            if idx < 1
                return result
            end

            nexc -= nvec[idx + 1] - 1
            nvec[idx + 1] = 0
            nvec[idx] += 1
            if nvec[idx] < n_max[idx]
                push!(result, copy(nvec))
            end
        end
    end
end

function _Importance(B::Vector{T}, bathPtr::AbstractVector, nvec::Vector{Int}) where T <: AbstractBath
    sum_γ = 0.0
    value = 1.0 + 0.0im
    
    for idx in findall(n_exc -> n_exc > 0, nvec)
        for _ in 1:(nvec[idx])
            k, ν = bathPtr[idx]
            γ = real(B[k][ν].γ)
            sum_γ += γ
            
            value *= ( B[k][ν].η / (γ * sum_γ) )
        end
    end 
    return abs(value)
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
        n_max = fill((tier + 1), Nterm)
    
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
            n_max = fill(1, Nterm)
        elseif tier >= 1
            n_max = fill(2, Nterm)
        end
    end

    # create idx2nvec and remove nvec when its value of importance is below threshold
    idx2nvec = _Idx2Nvec(n_max, tier)
    if threshold > 0.0
        drop_idx = Int[]
        for (idx, nvec) in enumerate(idx2nvec)            
            # only neglect the nvec where level ≥ 2
            level = sum(nvec)
            if level >= 2
                Ath = _Importance(B, bathPtr, nvec)
                if Ath < threshold
                    push!(drop_idx, idx)
                end
            end
        end
        deleteat!(idx2nvec, drop_idx)
    end

    # create lvl2idx and nvec2idx
    lvl2idx  = Dict{Int, Vector{Int}}()
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
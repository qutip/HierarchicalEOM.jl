abstract type AbstractHierarchyDict end

"""
    struct HierarchyDict <: AbstractHierarchyDict
An object which contains all dictionaries for pure (bosonic or fermionic) bath-ADOs hierarchy.

# Fields
- `idx2nvec` : Return the `Nvec` from a given index
- `nvec2idx` : Return the index from a given `Nvec`
- `lvl2idx` : Return the list of indices from a given level (excitation)
- `bathPtr` : Records the tuple ``(k, \\nu)`` for each position in `Nvec`, where ``k`` and ``\\nu`` represents the ``\\nu``-th exponential-expansion term of the ``k``-th bath.
"""
struct HierarchyDict <: AbstractHierarchyDict
    idx2nvec::Vector{Nvec}
    nvec2idx::Dict{Nvec, Int}
    lvl2idx::Dict{Int, Vector{Int}}
    bathPtr::Vector{Tuple}
end

"""
    struct MixHierarchyDict <: AbstractHierarchyDict
An object which contains all dictionaries for mixed (bosonic and fermionic) bath-ADOs hierarchy.

# Fields
- `idx2nvec` : Return the tuple `(Nvec_b, Nvec_f)` from a given index, where `b` represents boson and `f` represents fermion
- `nvec2idx` : Return the index from a given tuple `(Nvec_b, Nvec_f)`, where `b` represents boson and `f` represents fermion
- `Blvl2idx` : Return the list of indices from a given bosonic level (excitation)
- `Flvl2idx` : Return the list of indices from a given fermionic level (excitation)
- `bosonPtr` : Records the tuple ``(k, \\nu)`` for each position in `Nvec_b`, where ``k`` and ``\\nu`` represents the ``\\nu``-th exponential-expansion term of the ``k``-th bosonic bath.
- `fermionPtr` : Records the tuple ``(k, \\nu)`` for each position in `Nvec_f`, where ``k`` and ``\\nu`` represents the ``\\nu``-th exponential-expansion term of the ``k``-th fermionic bath.
"""
struct MixHierarchyDict <: AbstractHierarchyDict
    idx2nvec::Vector{Tuple{Nvec, Nvec}}
    nvec2idx::Dict{Tuple{Nvec, Nvec}, Int}
    Blvl2idx::Dict{Int, Vector{Int}}
    Flvl2idx::Dict{Int, Vector{Int}}
    bosonPtr::Vector{Tuple}
    fermionPtr::Vector{Tuple}
end

# generate index to n vector
function _Idx2Nvec(n_max::Vector{Int}, N_exc::Int)
    len = length(n_max)
    nvec = zeros(Int, len)
    result = [Nvec(nvec)]
    nexc = 0

    while true
        idx = len
        nvec[end] += 1
        nexc += 1
        if nvec[idx] < n_max[idx]
            push!(result, Nvec(nvec))
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
                push!(result, Nvec(nvec))
            end
        end
    end
end

function _Importance(B::Vector{T}, bathPtr::AbstractVector, nvec::Nvec) where T <: AbstractBath
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

# for pure hierarchy dictionary
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
            if nvec.level >= 2
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
    nvec2idx = Dict{Nvec, Int}()
    for level in 0:tier
        lvl2idx[level] = []
    end
    for (idx, nvec) in enumerate(idx2nvec)
        push!(lvl2idx[nvec.level], idx)
        nvec2idx[nvec] = idx
    end

    hierarchy = HierarchyDict(idx2nvec, nvec2idx, lvl2idx, bathPtr)
    return length(idx2nvec), baths, hierarchy
end

# for max hierarchy dictionary
function genBathHierarchy(bB::Vector{BosonBath}, fB::Vector{FermionBath}, tier_b::Int, tier_f::Int, dim::Int; threshold::Real=0.0)
    # deal with boson bath
    Nterm_b   = 0
    bosonPtr = Tuple[]
    baths_b = AbstractBosonBath[]
    for (k, b) in enumerate(bB)
        if b.dim != dim 
            error("The matrix size of the bosonic bath coupling operators are not consistent.")
        end
        push!(baths_b, b.bath...)
        for ν in 1:b.Nterm
            push!(bosonPtr, (k, ν))
        end
        Nterm_b += b.Nterm
    end
    n_max_b = fill((tier_b + 1), Nterm_b)
    idx2nvec_b = _Idx2Nvec(n_max_b, tier_b)

    # deal with fermion bath
    Nterm_f   = 0
    fermionPtr = Tuple[]
    baths_f = AbstractFermionBath[]
    for (k, b) in enumerate(fB)
        if b.dim != dim 
            error("The matrix size of the fermionic bath coupling operators are not consistent.")
        end
        push!(baths_f, b.bath...)
        for ν in 1:b.Nterm
            push!(fermionPtr, (k, ν))
        end
        Nterm_f += b.Nterm
    end
    if tier_f == 0
        n_max_f = fill(1, Nterm_f)
    elseif tier_f >= 1
        n_max_f = fill(2, Nterm_f)
    end
    idx2nvec_f = _Idx2Nvec(n_max_f, tier_f)

    # only store nvec tuple when its value of importance is above threshold
    idx2nvec = Tuple{Nvec, Nvec}[]
    if threshold > 0.0
        for nvec_b in idx2nvec_b
            for nvec_f in idx2nvec_f
                # only neglect the nvec tuple where level ≥ 2
                if (nvec_b.level >= 2) || (nvec_f.level >= 2)
                    Ath = _Importance(bB, bosonPtr, nvec_b) * _Importance(fB, fermionPtr, nvec_f)
                    if Ath >= threshold
                        push!(idx2nvec, (nvec_b, nvec_f))
                    end
                else
                    push!(idx2nvec, (nvec_b, nvec_f))
                end
            end
        end
    else
        for nvec_b in idx2nvec_b
            for nvec_f in idx2nvec_f
                push!(idx2nvec, (nvec_b, nvec_f))
            end
        end
    end

    # create lvl2idx and nvec2idx
    nvec2idx = Dict{Tuple{Nvec, Nvec}, Int}()
    blvl2idx = Dict{Int, Vector{Int}}()
    flvl2idx = Dict{Int, Vector{Int}}()
    for level in 0:tier_b
        blvl2idx[level] = []
    end
    for level in 0:tier_f
        flvl2idx[level] = []
    end
    for (idx, nvecTuple) in enumerate(idx2nvec)
        push!(blvl2idx[nvecTuple[1].level], idx)
        push!(flvl2idx[nvecTuple[2].level], idx)
        nvec2idx[nvecTuple] = idx
    end

    hierarchy = MixHierarchyDict(idx2nvec, nvec2idx, blvl2idx, flvl2idx, bosonPtr, fermionPtr)
    return length(idx2nvec), baths_b, baths_f, hierarchy
end
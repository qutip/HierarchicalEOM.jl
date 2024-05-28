abstract type AbstractHierarchyDict end

@doc raw"""
    struct HierarchyDict <: AbstractHierarchyDict
An object which contains all dictionaries for pure (bosonic or fermionic) bath-ADOs hierarchy.

# Fields
- `idx2nvec` : Return the `Nvec` from a given index of ADO
- `nvec2idx` : Return the index of ADO from a given `Nvec`
- `lvl2idx` : Return the list of ADO-indices from a given hierarchy level
- `bathPtr` : Records the tuple ``(\alpha, k)`` for each position in `Nvec`, where ``\alpha`` and ``k`` represents the ``k``-th exponential-expansion term of the ``\alpha``-th bath.
"""
struct HierarchyDict <: AbstractHierarchyDict
    idx2nvec::Vector{Nvec}
    nvec2idx::Dict{Nvec,Int}
    lvl2idx::Dict{Int,Vector{Int}}
    bathPtr::Vector{Tuple{Int,Int}}
end

@doc raw"""
    struct MixHierarchyDict <: AbstractHierarchyDict
An object which contains all dictionaries for mixed (bosonic and fermionic) bath-ADOs hierarchy.

# Fields
- `idx2nvec` : Return the tuple `(Nvec_b, Nvec_f)` from a given index of ADO, where `b` represents boson and `f` represents fermion
- `nvec2idx` : Return the index from a given tuple `(Nvec_b, Nvec_f)`, where `b` represents boson and `f` represents fermion
- `Blvl2idx` : Return the list of ADO-indices from a given bosonic-hierarchy level
- `Flvl2idx` : Return the list of ADO-indices from a given fermionic-hierarchy level
- `bosonPtr` : Records the tuple ``(\alpha, k)`` for each position in `Nvec_b`, where ``\alpha`` and ``k`` represents the ``k``-th exponential-expansion term of the ``\alpha``-th bosonic bath.
- `fermionPtr` : Records the tuple ``(\alpha, k)`` for each position in `Nvec_f`, where ``\alpha`` and ``k`` represents the ``k``-th exponential-expansion term of the ``\alpha``-th fermionic bath.
"""
struct MixHierarchyDict <: AbstractHierarchyDict
    idx2nvec::Vector{Tuple{Nvec,Nvec}}
    nvec2idx::Dict{Tuple{Nvec,Nvec},Int}
    Blvl2idx::Dict{Int,Vector{Int}}
    Flvl2idx::Dict{Int,Vector{Int}}
    bosonPtr::Vector{Tuple{Int,Int}}
    fermionPtr::Vector{Tuple{Int,Int}}
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

            nexc -= nvec[idx+1] - 1
            nvec[idx+1] = 0
            nvec[idx] += 1
            if nvec[idx] < n_max[idx]
                push!(result, Nvec(nvec))
            end
        end
    end
end

function _Importance(B::Vector{T}, bathPtr::AbstractVector, nvec::Nvec) where {T<:AbstractBath}
    sum_γ = 0.0
    value = 1.0 + 0.0im

    for idx in findall(n_exc -> n_exc > 0, nvec)
        for _ in 1:(nvec[idx])
            α, k = bathPtr[idx]
            γ = real(B[α][k].γ)
            sum_γ += γ

            value *= (B[α][k].η / (γ * sum_γ))
        end
    end
    return abs(value)
end

# for pure hierarchy dictionary
@noinline function genBathHierarchy(B::Vector{T}, tier::Int, dim::Int; threshold::Real = 0.0) where {T<:AbstractBath}
    Nterm = 0
    bathPtr = Tuple[]

    if T == BosonBath
        baths = AbstractBosonBath[]
        for (α, b) in enumerate(B)
            if b.dim != dim
                error("The matrix size of the bosonic bath coupling operators are not consistent.")
            end
            push!(baths, b.bath...)
            for k in 1:b.Nterm
                push!(bathPtr, (α, k))
            end
            Nterm += b.Nterm
        end
        n_max = fill((tier + 1), Nterm)

    elseif T == FermionBath
        baths = AbstractFermionBath[]
        for (α, b) in enumerate(B)
            if b.dim != dim
                error("The matrix size of the fermionic bath coupling operators are not consistent.")
            end
            push!(baths, b.bath...)
            for k in 1:b.Nterm
                push!(bathPtr, (α, k))
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
        splock = SpinLock()
        drop_idx = Int[]
        @threads for idx in eachindex(idx2nvec)
            nvec = idx2nvec[idx]

            # only neglect the nvec where level ≥ 2
            if nvec.level >= 2
                Ath = _Importance(B, bathPtr, nvec)
                if Ath < threshold
                    lock(splock)
                    try
                        push!(drop_idx, idx)
                    finally
                        unlock(splock)
                    end
                end
            end
        end
        sort!(drop_idx)
        deleteat!(idx2nvec, drop_idx)
    end

    # create lvl2idx and nvec2idx
    lvl2idx = Dict{Int,Vector{Int}}()
    nvec2idx = Dict{Nvec,Int}()
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

# for mixed hierarchy dictionary
@noinline function genBathHierarchy(
    bB::Vector{BosonBath},
    fB::Vector{FermionBath},
    tier_b::Int,
    tier_f::Int,
    dim::Int;
    threshold::Real = 0.0,
)
    # deal with boson bath
    Nterm_b = 0
    bosonPtr = Tuple[]
    baths_b = AbstractBosonBath[]
    for (α, b) in enumerate(bB)
        if b.dim != dim
            error("The matrix size of the bosonic bath coupling operators are not consistent.")
        end
        push!(baths_b, b.bath...)
        for k in 1:b.Nterm
            push!(bosonPtr, (α, k))
        end
        Nterm_b += b.Nterm
    end
    n_max_b = fill((tier_b + 1), Nterm_b)
    idx2nvec_b = _Idx2Nvec(n_max_b, tier_b)

    # deal with fermion bath
    Nterm_f = 0
    fermionPtr = Tuple[]
    baths_f = AbstractFermionBath[]
    for (α, b) in enumerate(fB)
        if b.dim != dim
            error("The matrix size of the fermionic bath coupling operators are not consistent.")
        end
        push!(baths_f, b.bath...)
        for k in 1:b.Nterm
            push!(fermionPtr, (α, k))
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
    idx2nvec = Tuple{Nvec,Nvec}[]
    if threshold > 0.0
        splock = SpinLock()
        @threads for nvec_b in idx2nvec_b
            for nvec_f in idx2nvec_f
                # only neglect the nvec tuple where level ≥ 2
                if (nvec_b.level >= 2) || (nvec_f.level >= 2)
                    Ath = _Importance(bB, bosonPtr, nvec_b) * _Importance(fB, fermionPtr, nvec_f)
                    if Ath >= threshold
                        lock(splock)
                        try
                            push!(idx2nvec, (nvec_b, nvec_f))
                        finally
                            unlock(splock)
                        end
                    end
                else
                    lock(splock)
                    try
                        push!(idx2nvec, (nvec_b, nvec_f))
                    finally
                        unlock(splock)
                    end
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
    nvec2idx = Dict{Tuple{Nvec,Nvec},Int}()
    blvl2idx = Dict{Int,Vector{Int}}()
    flvl2idx = Dict{Int,Vector{Int}}()
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

@doc raw"""
    getIndexEnsemble(nvec, bathPtr)
Search for all the multi-index ensemble ``(\alpha, k)`` where ``\alpha`` and ``k`` represents the ``k``-th exponential-expansion term in the ``\alpha``-th bath.

# Parameters
- `nvec::Nvec` : An object which records the repetition number of each multi-index ensembles in ADOs.
- `bathPtr::Vector{Tuple{Int, Int}}`: This can be obtained from [`HierarchyDict.bathPtr`](@ref HierarchyDict), [`MixHierarchyDict.bosonPtr`](@ref MixHierarchyDict), or [`MixHierarchyDict.fermionPtr`](@ref MixHierarchyDict).

# Returns
- `Vector{Tuple{Int, Int, Int}}`: a vector (list) of the tuples ``(\alpha, k, n)``.

# Example
Here is an example to use [`Bath`](@ref lib-Bath), [`Exponent`](@ref), [`HierarchyDict`](@ref), and `getIndexEnsemble` together:
```julia
L::M_Fermion;          # suppose this is a fermion type of HEOM Liouvillian superoperator matrix you create
HDict = L.hierarchy;   # the hierarchy dictionary
ados = SteadyState(L); # the stationary state (ADOs) for L 

# Let's consider all the ADOs for first level
idx_list = HDict.lvl2idx[1];

for idx in idx_list
    ρ1 = ados[idx]  # one of the 1st-level ADO
    nvec = HDict.idx2nvec[idx]  # the nvec corresponding to ρ1
    
    for (α, k, n) in getIndexEnsemble(nvec, HDict.bathPtr)
        α  # index of the bath
        k  # the index of the exponential-expansion term in α-th bath
        n  # the repetition number of the ensemble {α, k} in ADOs
        exponent = L.bath[α][k]  # the k-th exponential-expansion term in α-th bath

        # do some calculations you want
    end
end
```
"""
function getIndexEnsemble(nvec::Nvec, bathPtr::Vector{Tuple{Int,Int}})
    if length(nvec) != length(bathPtr)
        error("The given \"nvec\" and \"bathPtr\" are not consistent.")
    end

    result = Tuple{Int,Int,Int}[]
    for idx in nvec.data.nzind
        α, k = bathPtr[idx]
        push!(result, (α, k, nvec[idx]))
    end
    return result
end

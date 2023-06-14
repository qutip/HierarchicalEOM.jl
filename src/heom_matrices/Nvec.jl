@doc raw"""
    struct Nvec
An object which describes the repetition number of each multi-index ensembles in auxiliary density operators.

The `n_vector` (``\vec{n}``) denotes a set of integers:
```math
\{ n_{1,1}, ..., n_{\alpha, k}, ... \}
```
associated with the ``k``-th exponential-expansion term in the ``\alpha``-th bath.
If ``n_{\alpha, k} = 3`` means that the multi-index ensemble ``\{\alpha, k\}`` appears three times in the multi-index vector of ADOs (see the notations in our paper).

The hierarchy level (``L``) for an `n_vector` is given by ``L=\sum_{\alpha, k} n_{\alpha, k}``

# Fields
- `data` : the `n_vector`
- `level` : The level `L` for the `n_vector`

# Methods
One can obtain the repetition number for specific index (`idx`) by calling : `n_vector[idx]`.
To obtain the corresponding tuple ``(\alpha, k)`` for a given index `idx`, see `bathPtr` in [`HierarchyDict`](@ref) for more details.

`HEOM.jl` also supports the following calls (methods) :
```julia
length(n_vector);  # returns the length of `Nvec`
n_vector[1:idx];   # returns a vector which contains the excitation number of `n_vector` from index `1` to `idx`
n_vector[1:end];   # returns a vector which contains all the excitation number of `n_vector`
n_vector[:];       # returns a vector which contains all the excitation number of `n_vector`
from n in n_vector  # iteration
    # do something
end
```
"""
mutable struct Nvec
    data::SparseVector{Int, Int}
    level::Int
end

Nvec(V::Vector{Int}) = Nvec(sparsevec(V), sum(V))
Nvec(V::SparseVector{Int, Int}) = Nvec(copy(V), sum(V))

length(nvec::Nvec) = length(nvec.data)
lastindex(nvec::Nvec) = length(nvec)

getindex(nvec::Nvec, i::T) where {T <: Any} = nvec.data[i]

keys(nvec::Nvec) = keys(nvec.data)

function iterate(nvec::Nvec) 
    return nvec[1], 2
end
function iterate(nvec::Nvec, state::Int) 
    if state < length(nvec)
        return nvec[state], state + 1
    else
        return nvec[state], nothing
    end
end
iterate(nvec::Nvec, ::Nothing) = nothing

function show(io::IO, nvec::Nvec)
    print(io, "Nvec($(nvec[:]))")
end
function show(io::IO, m::MIME"text/plain", nvec::Nvec) show(io, nvec) end

hash(nvec::Nvec, h::UInt) = hash(nvec.data, h)
(==)(nvec1::Nvec, nvec2::Nvec) = hash(nvec1) == hash(nvec2)

copy(nvec::Nvec) = Nvec(copy(nvec.data), nvec.level)

function Nvec_plus!(nvec::Nvec, idx::Int) 
    nvec.data[idx] += 1
    nvec.level += 1
end
function Nvec_minus!(nvec::Nvec, idx::Int) 
    nvec.data[idx] -= 1
    nvec.level -= 1
end
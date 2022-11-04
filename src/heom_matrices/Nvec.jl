"""
    struct Nvec
An object which describes ``\\vec{n}``

The `n_vector` (``\\vec{n}``) denotes a set of integers:
```math
\\{ n_{11}, ..., n_{\\nu k}, ... \\}
```
where ``n_{\\nu k} \\geq 0`` is the excitation number associated with the ``k``-th exponential-expansion term in the ``\\nu``-th bath.

The hierarchy level (``L``) for an `n_vector` is given by ``L=\\sum_{\\nu, k} n_{\\nu k}``

# Fields
- `data` : the `n_vector`
- `level` : The level `L` for the `n_vector`

# Methods
One can obtain the excitation number for specific index (`idx`) by calling : `n_vector[idx]`.
To obtain the corresponding tuple ``(k, \\nu)`` for a given index `idx`, see `bathPtr` in [`HierarchyDict`](@ref) for more details.

`Heom.jl` also supports the following calls (methods) :
```julia
length(n_vector);  # returns the length of `Nvec`
n_vector[1:idx];   # returns a vector which contains the excitation number of `n_vector` from index `1` to `idx`
n_vector[1:end];   # returns a vector which contains all the excitation number of `n_vector`
n_vector[:];       # returns a vector which contains all the excitation number of `n_vector`
from n in n_vector  # iteration
    # do something
end
"""
mutable struct Nvec
    data::Vector{Int}
    level::Int
end

Nvec(V::Vector{Int}) = Nvec(copy(V), sum(V))

length(nvec::Nvec) = length(nvec.data)
lastindex(nvec::Nvec) = length(nvec)

getindex(nvec::Nvec, i::T) where {T <: Any} = nvec.data[i]

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
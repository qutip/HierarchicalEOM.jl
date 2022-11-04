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
- `value` : In stead of storing the exact vector, we transform the `nvec` into a `value`
- `level` : The level `L` for the `n_vector`
- `base` : equals to `maximum_tier + 1`, which is used to transfrom between `n_vector` and `value`
- `digit` : equals to the length of `n_vector`, which is used to transfrom between `n_vector` and `value`

# Methods
One can obtain the excitation number for specific index (`idx`) by calling : `n_vector[idx]`.
To obtain the corresponding tuple ``(k, \\nu)`` for a given index `idx`, see `bathPtr` in [`HierarchyDict`](@ref) for more details.

`Heom.jl` also supports the following calls (methods) :
```julia
length(n_vector);  # returns the length of `Nvec`
ados[1:idx];   # returns a vector which contains the excitation number of `n_vector` from index `1` to `idx`
ados[1:end];   # returns a vector which contains all the excitation number of `n_vector`
ados[:];       # returns a vector which contains all the excitation number of `n_vector`
from n in n_vector  # iteration
    # do something
end
"""
mutable struct Nvec
    value::Int
    level::Int
    const base::Int
    const digit::Int
end

function Nvec(V::Vector{Int}, base::Int)
    value = 0
    level = 0
    L = length(V)
    for i in findall(ai -> ai > 0, V)
        value += V[i] * base ^ (L - i)
        level += V[i]
    end
        
    return Nvec(value, level, base, L)
end

function checkbounds(nvec::Nvec, i::Int)
    if (i > nvec.digit) || (i < 1)
        error("Attempt to access $(nvec.digit)-element Nvec at index [$(i)]")
    end
end

length(nvec::Nvec) = nvec.digit
lastindex(nvec::Nvec) = length(nvec)

function getindex(nvec::Nvec, i::Int)
    checkbounds(nvec, i)
    return (nvec.value ÷ nvec.base ^ (nvec.digit - i)) % nvec.base
end
function getindex(nvec::Nvec, r::UnitRange{Int})
    if length(r) > 0
        checkbounds(nvec, r[1])
        checkbounds(nvec, r[end])

        result = zeros(Int, length(r))
        
        basis = nvec.base ^ (nvec.digit - r[1])
        result[1] = (nvec.value ÷ basis) % nvec.base
        value = nvec.value % basis
        for i in 2:length(r)
            if value == 0
                break
            else
                basis ÷= nvec.base
                result[i] = value ÷ basis
                value %= basis
            end
        end
        return result
    else
        return Int[]
    end
end
getindex(nvec::Nvec, ::Colon) = getindex(nvec, 1:lastindex(nvec))

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

hash(nvec::Nvec, h::UInt) = hash((nvec.value, nvec.base, nvec.digit), h)
(==)(nvec1::Nvec, nvec2::Nvec) = hash(nvec1) == hash(nvec2)

copy(nvec::Nvec) = Nvec(nvec.value, nvec.level, nvec.base, nvec.digit)

function prev_grad(nvec::Nvec, i::Int)
    basis = nvec.base ^ (nvec.digit - i)
    if (nvec.value ÷ basis) % nvec.base > 0  # nvec[i] > 0
        return basis
    else
        return 0
    end
end

function next_grad(nvec::Nvec, i::Int, tier::Int)
    if nvec.level < tier
        return nvec.base ^ (nvec.digit - i)
    else
        0
    end
end

function Nvec_plus!(nvec::Nvec, value::Int) 
    nvec.value += value
    nvec.level += 1
end
function Nvec_minus!(nvec::Nvec, value::Int) 
    nvec.value -= value
    nvec.level -= 1
end
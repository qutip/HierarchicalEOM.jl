"""
    mutable struct ADOs
The Auxiliary Density Operators for Heom model.

# Fields
- `data` : the vectorized auxiliary density operators
- `dim` : the dimension of the system
- `N` : the number of auxiliary density operators

# Methods
One can obtain the density matrix for specific index (`idx`) by calling : `ados[idx]`.
`Heom.jl` also supports the following calls (methods) :
```julia
length(ados);  # returns the total number of `ADOs`
ados[1:idx];   # returns a vector which contains the `ADO` (in matrix form) from index `1` to `idx`
ados[1:end];   # returns a vector which contains all the `ADO` (in matrix form)
ados[:];       # returns a vector which contains all the `ADO` (in matrix form)
from rho in ados  # iteration
    # do something
end
```
"""
mutable struct ADOs 
    data::SparseVector{ComplexF64, Int64}
    const dim::Int
    const N::Int
end

"""
    ADOs(V, N)
Gernerate the object of auxiliary density operators for Heom model.

# Parameters
- `V::AbstractVector` : the vectorized auxiliary density operators
- `N::Int` : the number of auxiliary density operators.
"""
function ADOs(V::AbstractVector, N::Int)
    # check the dimension of V
    d,  = size(V)
    dim = √(d / N)
    if isinteger(dim)
        return ADOs(sparsevec(V), Int(dim), N)
    else
        error("The dimension of vector is not consistent with N.")
    end    
end

function checkbounds(A::ADOs, i::Int)
    if (i > A.N) || (i < 1)
        error("Attempt to access $(A.N)-element ADOs at index [$(i)]")
    end
end

"""
    length(A::ADOs)
Returns the total number of the Auxiliary Density Operators (ADOs)
"""
length(A::ADOs) = A.N

lastindex(A::ADOs) = length(A)

function getindex(A::ADOs, i::Int)
    checkbounds(A, i)

    sup_dim = A.dim ^ 2
    back    = sup_dim * i
    return sparse(reshape(A.data[(back - sup_dim + 1):back], A.dim, A.dim))
end

function getindex(A::ADOs, r::UnitRange{Int})
    checkbounds(A, r[1])
    checkbounds(A, r[end])

    result = []
    sup_dim = A.dim ^ 2
    for i in r
        back = sup_dim * i
        push!(result, 
            sparse(reshape(A.data[(back - sup_dim + 1):back], A.dim, A.dim))
        )
    end
    return result
end
getindex(A::ADOs, ::Colon) = getindex(A, 1:lastindex(A))

function iterate(A::ADOs) 
    return A[1], 2
end
function iterate(A::ADOs, state::Int) 
    if state < length(A)
        return A[state], state + 1
    else
        return A[state], nothing
    end
end
iterate(A::ADOs, ::Nothing) = nothing

function show(io::IO, A::ADOs)
    print(io, 
        "Auxiliary Density Operators with (system) dim = $(A.dim), N = $(A.N)\n"
    )
end
function show(io::IO, m::MIME"text/plain", A::ADOs) show(io, A) end

"""
    getRho(ados)
Return the density matrix of the reduced state (system) from a given auxiliary density operators

# Parameters
- `ados::ADOs` : the auxiliary density operators for Heom model

# Returns
- `ρ` : The density matrix of the reduced state
"""
function getRho(ados::ADOs)
    return sparse(reshape(ados.data[1:((ados.dim) ^ 2)], ados.dim, ados.dim))
end

"""
    getADO(ados, idx)
Return the auxiliary density operator with a specific index from auxiliary density operators

This function equals to calling : `ados[idx]`.

# Parameters
- `ados::ADOs` : the auxiliary density operators for Heom model
- `idx::Int` : the index of the auxiliary density operator

# Returns
- `ρ_idx` : The auxiliary density operator
"""
getADO(ados::ADOs, idx::Int) = ados[idx]
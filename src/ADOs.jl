"""
    mutable struct ADOs
The Auxiliary Density Operators for Heom model.

# Fields
- `data` : the vectorized auxiliary density operators
- `dim` : the dimension of the system
- `Nb` : the number of bosonic states
- `Nf` : the number of fermionic states

# Methods
For pure bosonic (or fermionic) type bath `ADOs`, 
one can obtain the density matrix for specific index (`idx`) by calling : `ados[idx]`.
`Heom.jl` also supports the following calls (methods) :
```@example
length(ados);  # returns the total number of `ADOs`
ados[1:idx];   # returns a vector which contains the `ADO` (in matrix form) from index `1` to `idx`
ados[1:end];   # returns a vector which contains all the `ADO` (in matrix form)
ados[:];       # returns a vector which contains all the `ADO` (in matrix form)
from rho in ados  # iteration
    # do something
end
```

For mixed (bosonic and fermionic) type bath `ADOs`, 
one needs two indices (`idx_b` and `idx_f`), and thus, `ados[idx_b, idx_f]`.
Note that the first index specifies the bosonic bath index while the other one specifies the fermionic bath.
`Heom.jl` also supports the following calls (methods) :
```@example
length(ados);       # returns the total number of `ADOs`
ados[idx_b, 1:end]; # returns a vector which contains all the fermionic `ADO` (in matrix form) where bosonic index is `idx_b`
ados[1:end, idx_f]; # returns a vector which contains all the bosonic `ADO` (in matrix form) where bosonic index is `idx_f`
ados[idx_b, :];     # returns a vector which contains all the fermionic `ADO` (in matrix form) where bosonic index is `idx_b`
ados[:, idx_f];     # returns a vector which contains all the bosonic `ADO` (in matrix form) where bosonic index is `idx_f`
```
But, currently, we don't support `iterate()` for mixed bath ADOs.
"""
mutable struct ADOs 
    data::SparseVector{ComplexF64, Int64}
    const dim::Int
    const Nb::Int
    const Nf::Int
end

"""
    ADOs(V, Nb=0, Nf=0)
Gernerate the object of auxiliary density operators for Heom model.

# Parameters
- `V::AbstractVector` : the vectorized auxiliary density operators
- `Nb::Int` : the number of bosonic states. Defaults to `0`.
- `Nf::Int` : the number of fermionic states Defaults to `0`.
"""
function ADOs(        
        V::AbstractVector,
        Nb::Int=0,
        Nf::Int=0
    )

    # get total number of Env. states
    local Ntot :: Int
    if Nb == 0     # fermionic bath
        Ntot = Nf

    elseif Nf == 0 # bosonic bath
        Ntot = Nb

    else               # mixed bath
        Ntot = Nb * Nf
    end

    # check the dimension of V
    d,  = size(V)
    dim = √(d / Ntot)
    if isinteger(dim)
        return ADOs(sparsevec(V), Int(dim), Nb, Nf)
    else
        error("The dimension of vector is not consistent with Nb and Nf.")
    end    
end

function checkbounds(A::ADOs, i::Int, tag::String)
    # boson case
    if tag == "b"
        if (i > A.Nb) || (i < 1)
            error("BoundsError: attempt to access $(A.Nb)-element (bosonic) ADOs at index [$(i)]")
        end

    # fermion case
    elseif tag == "f"
        if (i > A.Nf) || (i < 1)
            error("BoundsError: attempt to access $(A.Nf)-element (fermionic) ADOs at index [$(i)]")
        end
    end
end

"""
    length(A::ADOs)
Returns the total number of the Auxiliary Density Operators (ADOs)
"""
function length(A::ADOs)
    if A.Nb == 0
        return A.Nf
    elseif A.Nf == 0
        return A.Nb
    else
        return A.Nb * A.Nf
    end
end

lastindex(A::ADOs) = length(A)
function lastindex(A::ADOs, d::Int)
    if d == 1  # Boson
        return A.Nb
    elseif d == 2  # fermion
        return A.Nf
    end
end

function getindex(A::ADOs, i::Int)
    if A.Nb == 0 
        checkbounds(A, i, "f")
    elseif A.Nf == 0 
        checkbounds(A, i, "b")
    else
        error("The ADOs is from mixed (bosonic and fermionic) bath, use \"ados[idx_b, idx_f]\" or function \"getADO(ados, idx_b, idx_f)\" instead.")
    end

    sup_dim = A.dim ^ 2
    back    = sup_dim * i
    return sparse(reshape(A.data[(back - sup_dim + 1):back], A.dim, A.dim))
end

function getindex(A::ADOs, r::UnitRange{Int})
    if A.Nb == 0 
        checkbounds(A, r[1],   "f")
        checkbounds(A, r[end], "f")
    elseif A.Nf == 0 
        checkbounds(A, r[1],   "b")
        checkbounds(A, r[end], "b")
    else
        error("The ADOs is from mixed (bosonic and fermionic) bath, use \"ados[idx_b, idx_f]\" instead.")
    end

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

function getindex(A::ADOs, i::Int, j::Int)
    if A.Nb == 0
        error("The given ADOs is from pure fermionic bath, use \"ados[idx]\" or function \"getADO(ados, idx)\" instead.")
    
    elseif A.Nf == 0
        error("The given ADOs is from pure bosonic bath, use \"ados[idx]\" or function \"getADO(ados, idx)\" instead.")
    
    else
        checkbounds(A, i, "b")
        checkbounds(A, j, "f")
        
        idx     = (i - 1) * A.Nf + j
        sup_dim = A.dim ^ 2
        back    = sup_dim * idx
        return sparse(reshape(A.data[(back - sup_dim + 1):back], A.dim, A.dim))
    end
end

function getindex(A::ADOs, i::Int, r::UnitRange{Int})
    if A.Nb == 0
        error("The given ADOs is from pure fermionic bath, use \"ados[idx]\" instead.")
    
    elseif A.Nf == 0
        error("The given ADOs is from pure bosonic bath, use \"ados[idx]\" instead.")
    
    else
        checkbounds(A, i, "b")
        checkbounds(A, r[1],   "f")
        checkbounds(A, r[end], "f")
        
        result = []
        sup_dim = A.dim ^ 2
        idx_b = (i - 1) * A.Nf
        for j in r
            idx  = idx_b + j
            back = sup_dim * idx
            push!(result, 
                sparse(reshape(A.data[(back - sup_dim + 1):back], A.dim, A.dim))
            )
        end
        return result
    end
end

function getindex(A::ADOs, r::UnitRange{Int}, j::Int)
    if A.Nb == 0
        error("The given ADOs is from pure fermionic bath, use \"ados[idx]\" instead.")
    
    elseif A.Nf == 0
        error("The given ADOs is from pure bosonic bath, use \"ados[idx]\" instead.")
    
    else
        checkbounds(A, r[1],   "b")
        checkbounds(A, r[end], "b")
        checkbounds(A, j, "f")
    
        result = []
        sup_dim = A.dim ^ 2
        for i in r
            idx  = (i - 1) * A.Nf + j
            back = sup_dim * idx
            push!(result, 
                sparse(reshape(A.data[(back - sup_dim + 1):back], A.dim, A.dim))
            )
        end
        return result
    end
end
getindex(A::ADOs, ::Colon, j::Int) = getindex(A, 1:lastindex(A, 1), j)
getindex(A::ADOs, i::Int, ::Colon) = getindex(A, i, 1:lastindex(A, 2))

function iterate(A::ADOs) 
    if (A.Nb == 0) || (A.Nf == 0)
        return A[1], 2
    else
        error("ADOs doesn't support \"iterate\" for mixed (bosonic and fermionic) bath yet.")
    end
end
function iterate(A::ADOs, state) 
    if state < length(A)
        return A[state], state + 1
    else
        return A[state], nothing
    end
end
iterate(A::ADOs, ::Nothing) = nothing

function show(io::IO, A::ADOs)
    print(io, 
        "Auxiliary Density Operators with (system) dim = $(A.dim), Nb = $(A.Nb), Nf = $(A.Nf)\n"
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

"""
    getADO(ados, idx_b, idx_f)
Return the auxiliary density operator with specific indices *[only for mixtured (bosonic and fermionic) bath]*

This function equals to calling : `ados[idx_b, idx_f]`.

# Parameters
- `ados::ADOs` : the auxiliary density operators for Heom model
- `idx_b::Int` : the bosonic-state index of the auxiliary density operator.
- `idx_f::Int` : the fermionic-state index of the auxiliary density operator.

# Returns
- `ρ_idx` : The auxiliary density operator
"""
getADO(ados::ADOs, idx_b::Int, idx_f::Int) = ados[idx_b, idx_f]
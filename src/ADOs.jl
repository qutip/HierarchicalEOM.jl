"""
    mutable struct ADOs
The Auxiliary Density Operators for Heom model.

# Fields
- `data` : the vectorized auxiliary density operators
- `dim` : the dimension of the system
- `Nb` : the number of bosonic states
- `Nf` : the number of fermionic states
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

"""
    size(A::ADOs)
Returns the size of the vetorized Auxiliary Density Operators (ADOs)
"""
size(A::ADOs) = size(A.data)

function show(io::IO, V::ADOs)
    print(io, 
        "Auxiliary Density Operators with (system) dim = $(V.dim)\n",
        "bosonic-state   number Nb = $(V.Nb)\n",
        "fermionic-state number Nf = $(V.Nf)\n",
        "data =\n"
    )
    show(io, MIME("text/plain"), V.data)
end

function show(io::IO, m::MIME"text/plain", V::ADOs) show(io, V) end

"""
    getRho(ados)
Return the density matrix of the reduced state (system) from a given auxiliary density operators

# Parameters
- `ados::ADOs` : the auxiliary density operators for Heom model

# Returns
- `ρ` : The density matrix of the reduced state
"""
function getRho(ados::ADOs)
    return Matrix(reshape(ados.data[1:((ados.dim) ^ 2)], ados.dim, ados.dim))
end

"""
    getADO(ados, idx)
Return the auxiliary density operator with a specific index from auxiliary density operators

# Parameters
- `ados::ADOs` : the auxiliary density operators for Heom model
- `idx::Int` : the index of the auxiliary density operator

# Returns
- `ρ_idx` : The auxiliary density operator
"""
function getADO(ados::ADOs, idx::Int)
    if (ados.Nb != 0) && (ados.Nf != 0)
        error("The given ADOs is from mixed (bosonic and fermionic) bath, use function \"getADO(ados, idx_b, idx_f)\" instead.")
    end

    if idx < 1
        error("idx must be greater than 0")
    end

    if idx == 1
        return getRho(ados)
    end

    # check whether the index is for fermion or boson or mixed env.
    if ((ados.Nb == 0) && (idx <= ados.Nf)) || ((ados.Nf == 0) && (idx <= ados.Nb)) || (idx <= (ados.Nb * ados.Nf))
        sup_dim = ados.dim ^ 2
        back    = sup_dim * idx
        return Matrix(reshape(ados.data[(back - sup_dim + 1):back], ados.dim, ados.dim))
    
    else
        error("idx is not available.")
    end
end

"""
    getADO(ados, idx_b, idx_f)
Return the auxiliary density operator with specific indices *[only for mixtured (bosonic and fermionic) bath]*

# Parameters
- `ados::ADOs` : the auxiliary density operators for Heom model
- `idx_b::Int` : the bosonic-state index of the auxiliary density operator.
- `idx_f::Int` : the fermionic-state index of the auxiliary density operator.

# Returns
- `ρ_idx` : The auxiliary density operator
"""
function getADO(ados::ADOs, idx_b::Int, idx_f::Int)
    if ados.Nb == 0
        error("The given ADOs is from pure fermionic bath, use function \"getADO(ados, idx)\" instead.")
    
    elseif ados.Nf == 0
        error("The given ADOs is from pure bosonic bath, use function \"getADO(ados, idx)\" instead.")
    
    else
        if (idx_b == 1) && (idx_f == 1)
            return getRho(ados)

        elseif (idx_b <= ados.Nb) && (idx_f <= ados.Nf)
            idx     = (idx_b - 1) * ados.Nf + idx_f
            sup_dim = ados.dim ^ 2
            back    = sup_dim * idx
            return Matrix(reshape(ados.data[(back - sup_dim + 1):back], ados.dim, ados.dim))

        else
            error("idx_b and idx_f are not available.")
        end 
    end
end
"""
# `ADOs`
The Auxiliary Density Operators for Heom model.

## Fields
- `data::Vector{ComplexF64}` : the sparse vector
- `dim::Int` : the dimension of the system
- `Nb::Int` : the number of bosonic states
- `Nf::Int` : the number of fermionic states

## Constructor
`ADOs(V, dim, Nb, Nf)`

- `data::Vector{ComplexF64}` : the sparse vector
- `dim::Int` : the dimension of the system
- `Nb::Int` : the number of bosonic states. Defaults to 0.
- `Nf::Int` : the number of fermionic states Defaults to 0.
"""
mutable struct ADOs 
    data::Vector{ComplexF64}
    const dim::Int
    const Nb::Int
    const Nf::Int

    function ADOs(        
            V::AbstractVector,
            dim::Int,
            Nb::Int=0,
            Nf::Int=0
        )

        if dim < 2
            error("dim must be greater than 1")
        end

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
        d, = size(V)
        if d != (dim ^ 2) * Ntot
            error("The dimension is not consistent, the size of \"V\" should be equal to \"(dim ^ 2) * Nb * Nf\".")
        end

        # make eltype of V: ComplexF64
        if eltype(V) != ComplexF64
            V = convert.(ComplexF64, V)
        end
        return new(Vector(V), dim, Nb, Nf)
    end
end

size(A::ADOs) = size(A.data)

"""
# `getRho(ados)`
Return the density matrix of the reduced state (system) from a given auxiliary density operators

## Parameters
- `ados::ADOs` : the auxiliary density operators for Heom model

## Returns
- `ρ` : The density matrix of the reduced state
"""
function getRho(ados::ADOs)
    return sparse(reshape(ados.data[1:((ados.dim) ^ 2)], ados.dim, ados.dim))
end

"""
# `getADO(ados, idx)`
Return the auxiliary density operator with a specific index from auxiliary density operators

## Parameters
- `ados::ADOs` : the auxiliary density operators for Heom model
- `idx::Int`   : the index of the auxiliary density operator

## Returns
- `ρ_idx` : The auxiliary density operator
"""
function getADO(ados::ADOs, idx::Int)
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
        return sparse(reshape(ados.data[(back - sup_dim + 1):back], ados.dim, ados.dim))
    
    else
        error("idx is not available.")
    end
end

"""
# `getADO(ados, idx_b, idx_f)`
Return the auxiliary density operator with specific indices *[only for mixtured (bosonic and fermionic) bath]*

## Parameters
- `ados::ADOs` : the auxiliary density operators for Heom model
- `idx_b::Int`   : the bosonic-state index of the auxiliary density operator. Defaults to 1
- `idx_f::Int`   : the fermionic-state index of the auxiliary density operator. Defaults to 1

## Returns
- `ρ_idx` : The auxiliary density operator
"""
function getADO(ados::ADOs, idx_b::Int=1, idx_f::Int=1)
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
            return sparse(reshape(ados.data[(back - sup_dim + 1):back], ados.dim, ados.dim))

        else
            error("idx_b and idx_f are not available.")
        end 
    end
end
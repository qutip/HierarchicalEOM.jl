@doc raw"""
    struct M_S <: AbstractHEOMMatrix
HEOM Liouvillian superoperator matrix with cutoff level of the hierarchy equals to `0`.  
This corresponds to the standard Schrodinger (Liouville-von Neumann) equation, namely
```math
M[\cdot]=-i \left[H_{sys}, \cdot \right]_-,
```
where ``[\cdot, \cdot]_-`` stands for commutator.

# Fields
- `data` : the sparse matrix of HEOM liouvillian superoperator
- `tier` : the tier (cutoff level) for the hierarchy, which equals to `0` in this case
- `dim` : the dimension of system
- `N` : the number of total ADOs, which equals to `1` (only the reduced density operator) in this case
- `sup_dim` : the dimension of system superoperator
- `parity` : the parity label of the fermionic system (usually `:even`, only set as `:odd` for calculating spectrum of fermionic system)
"""
struct M_S <: AbstractHEOMMatrix
    data::SparseMatrixCSC{ComplexF64, Int64}
    tier::Int
    dim::Int
    N::Int
    sup_dim::Int
    parity::Symbol
end

@doc raw"""
    M_S(Hsys, parity=:even; verbose=true)
Generate HEOM Liouvillian superoperator matrix with cutoff level of the hierarchy equals to `0`.  
This corresponds to the standard Schrodinger (Liouville-von Neumann) equation, namely
```math
M[\cdot]=-i \left[H_{sys}, \cdot \right]_-,
```
where ``[\cdot, \cdot]_-`` stands for commutator.

# Parameters
- `Hsys` : The time-independent system Hamiltonian
- `parity::Symbol` : the parity label of the fermionic system (only set as `:odd` for calculating spectrum of fermionic system). Defaults to `:even`.
- `verbose::Bool` : To display verbose output during the process or not. Defaults to `true`.

Note that the parity only need to be set as `:odd` when the system contains fermionic systems and you need to calculate the spectrum (density of states) of it.
"""
@noinline function M_S(        
        Hsys,
        parity::Symbol=:even;
        verbose::Bool=true
    )

    # check parity
    if (parity != :even) && (parity != :odd)
        error("The parity symbol of density matrix should be either \":even\" or \":odd\".")
    end

    # check for system dimension
    if !isValidMatrixType(Hsys)
        error("Invalid matrix \"Hsys\" (system Hamiltonian).")
    end
    Nsys,   = size(Hsys)
    sup_dim = Nsys ^ 2

    # the liouvillian operator for free Hamiltonian
    if verbose
        println("Constructing Liouville-von Neumann superoperator...")
        flush(stdout)
    end
    Lsys = -1im * (spre(Hsys) - spost(Hsys))
    if verbose
        println("[DONE]")
        flush(stdout)
    end
    return M_S(Lsys, 0, Nsys, 1, sup_dim, parity)
end
@doc raw"""
    struct M_S{T} <: AbstractHEOMLSMatrix
HEOM Liouvillian superoperator matrix with cutoff level of the hierarchy equals to `0`.  
This corresponds to the standard Schrodinger (Liouville-von Neumann) equation, namely
```math
M[\cdot]=-i \left[H_{sys}, \cdot \right]_-,
```
where ``[\cdot, \cdot]_-`` stands for commutator.

# Fields
- `data::T` : the sparse matrix of HEOM Liouvillian superoperator
- `tier` : the tier (cutoff level) for the hierarchy, which equals to `0` in this case
- `dim` : the dimension of system
- `N` : the number of total ADOs, which equals to `1` (only the reduced density operator) in this case
- `sup_dim` : the dimension of system superoperator
- `parity` : the parity label of the operator which HEOMLS is acting on (usually `EVEN`, only set as `ODD` for calculating spectrum of fermionic system).
"""
struct M_S{T} <: AbstractHEOMLSMatrix
    data::T
    tier::Int
    dim::Int
    N::Int
    sup_dim::Int
    parity::AbstractParity
end

@doc raw"""
    M_S(Hsys, parity=EVEN; verbose=true)
Generate HEOM Liouvillian superoperator matrix with cutoff level of the hierarchy equals to `0`.  
This corresponds to the standard Schrodinger (Liouville-von Neumann) equation, namely
```math
M[\cdot]=-i \left[H_{sys}, \cdot \right]_-,
```
where ``[\cdot, \cdot]_-`` stands for commutator.

# Parameters
- `Hsys` : The time-independent system Hamiltonian
- `parity::AbstractParity` : the parity label of the operator which HEOMLS is acting on (usually `EVEN`, only set as `ODD` for calculating spectrum of fermionic system).
- `verbose::Bool` : To display verbose output during the process or not. Defaults to `true`.

Note that the parity only need to be set as `ODD` when the system contains fermionic systems and you need to calculate the spectrum (density of states) of it.
"""
@noinline function M_S(        
        Hsys,
        parity::AbstractParity=EVEN;
        verbose::Bool=true
    )

    # check for system dimension
    _Hsys = HandleMatrixType(Hsys, 0, "Hsys (system Hamiltonian)")
    Nsys    = size(_Hsys, 1)
    sup_dim = Nsys ^ 2

    # the Liouvillian operator for free Hamiltonian
    if verbose
        println("Constructing Liouville-von Neumann superoperator...")
        flush(stdout)
    end
    Lsys = minus_i_L_op(_Hsys)
    if verbose
        println("[DONE]")
        flush(stdout)
    end
    return M_S{SparseMatrixCSC{ComplexF64, Int64}}(Lsys, 0, Nsys, 1, sup_dim, parity)
end
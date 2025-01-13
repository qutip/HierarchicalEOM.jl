export M_Fermion

@doc raw"""
    struct M_Fermion <: AbstractHEOMLSMatrix
HEOM Liouvillian superoperator matrix for fermionic bath

# Fields
- `data<:AbstractSciMLOperator` : the matrix of HEOM Liouvillian superoperator
- `tier` : the tier (cutoff level) for the fermionic hierarchy
- `dimensions` : the dimension list of the coupling operator (should be equal to the system dimensions).
- `N` : the number of total ADOs
- `sup_dim` : the dimension of system superoperator
- `parity` : the parity label of the operator which HEOMLS is acting on (usually `EVEN`, only set as `ODD` for calculating spectrum of fermionic system).
- `bath::Vector{FermionBath}` : the vector which stores all `FermionBath` objects
- `hierarchy::HierarchyDict`: the object which contains all dictionaries for fermion-bath-ADOs hierarchy.

!!! note "`dims` property"
    For a given `M::M_Fermion`, `M.dims` or `getproperty(M, :dims)` returns its `dimensions` in the type of integer-vector.
"""
struct M_Fermion{T<:AbstractSciMLOperator} <: AbstractHEOMLSMatrix{T}
    data::T
    tier::Int
    dimensions::Dimensions
    N::Int
    sup_dim::Int
    parity::AbstractParity
    bath::Vector{FermionBath}
    hierarchy::HierarchyDict
end

function M_Fermion(
    Hsys::QuantumObject,
    tier::Int,
    Bath::FermionBath,
    parity::AbstractParity = EVEN;
    threshold::Real = 0.0,
    verbose::Bool = true,
)
    return M_Fermion(Hsys, tier, [Bath], parity, threshold = threshold, verbose = verbose)
end

@doc raw"""
    M_Fermion(Hsys, tier, Bath, parity=EVEN; threshold=0.0, verbose=true)
Generate the fermion-type HEOM Liouvillian superoperator matrix

# Parameters
- `Hsys` : The time-independent system Hamiltonian or Liouvillian
- `tier::Int` : the tier (cutoff level) for the fermionic bath
- `Bath::Vector{FermionBath}` : objects for different fermionic baths
- `parity::AbstractParity` : the parity label of the operator which HEOMLS is acting on (usually `EVEN`, only set as `ODD` for calculating spectrum of fermionic system).
- `threshold::Real` : The threshold of the importance value (see Ref. [1]). Defaults to `0.0`.
- `verbose::Bool` : To display verbose output and progress bar during the process or not. Defaults to `true`.

[1] [Phys. Rev. B 88, 235426 (2013)](https://doi.org/10.1103/PhysRevB.88.235426)
"""
@noinline function M_Fermion(
    Hsys::QuantumObject,
    tier::Int,
    Bath::Vector{FermionBath},
    parity::AbstractParity = EVEN;
    threshold::Real = 0.0,
    verbose::Bool = true,
)

    # check for system dimension
    _Hsys = HandleMatrixType(Hsys, "Hsys (system Hamiltonian or Liouvillian)")
    sup_dim = prod(_Hsys.dimensions)^2
    I_sup = sparse(one(ComplexF64) * I, sup_dim, sup_dim)

    # the Liouvillian operator for free Hamiltonian term
    Lsys = minus_i_L_op(_Hsys)

    # fermionic bath
    if verbose && (threshold > 0.0)
        print("Checking the importance value for each ADOs...")
        flush(stdout)
    end
    Nado, baths, hierarchy = genBathHierarchy(Bath, tier, _Hsys.dimensions, threshold = threshold)
    idx2nvec = hierarchy.idx2nvec
    nvec2idx = hierarchy.nvec2idx
    if verbose && (threshold > 0.0)
        println("[DONE]")
        flush(stdout)
    end

    # start to construct the matrix
    Nthread = nthreads()
    L_row = [Int[] for _ in 1:Nthread]
    L_col = [Int[] for _ in 1:Nthread]
    L_val = [ComplexF64[] for _ in 1:Nthread]

    if verbose
        println("Preparing block matrices for HEOM Liouvillian superoperator (using $(Nthread) threads)...")
        flush(stdout)
        prog = ProgressBar(Nado)
    end
    @threads for idx in 1:Nado
        tID = threadid()

        # fermion (current level) superoperator
        nvec = idx2nvec[idx]
        if nvec.level >= 1
            sum_γ = bath_sum_γ(nvec, baths)
            op = Lsys - sum_γ * I_sup
        else
            op = Lsys
        end
        add_operator!(op, L_row[tID], L_col[tID], L_val[tID], Nado, idx, idx)

        # connect to fermionic (n+1)th- & (n-1)th- level superoperator
        mode = 0
        nvec_neigh = copy(nvec)
        for fB in baths
            for k in 1:fB.Nterm
                mode += 1
                n_k = nvec[mode]

                # connect to fermionic (n-1)th-level superoperator
                if n_k > 0
                    Nvec_minus!(nvec_neigh, mode)
                    if (threshold == 0.0) || haskey(nvec2idx, nvec_neigh)
                        idx_neigh = nvec2idx[nvec_neigh]
                        op = minus_i_C_op(fB, k, nvec.level, sum(nvec_neigh[1:(mode-1)]), parity)
                        add_operator!(op, L_row[tID], L_col[tID], L_val[tID], Nado, idx, idx_neigh)
                    end
                    Nvec_plus!(nvec_neigh, mode)

                    # connect to fermionic (n+1)th-level superoperator
                elseif nvec.level < tier
                    Nvec_plus!(nvec_neigh, mode)
                    if (threshold == 0.0) || haskey(nvec2idx, nvec_neigh)
                        idx_neigh = nvec2idx[nvec_neigh]
                        op = minus_i_A_op(fB, nvec.level, sum(nvec_neigh[1:(mode-1)]), parity)
                        add_operator!(op, L_row[tID], L_col[tID], L_val[tID], Nado, idx, idx_neigh)
                    end
                    Nvec_minus!(nvec_neigh, mode)
                end
            end
        end
        if verbose
            next!(prog) # trigger a progress bar update
        end
    end
    if verbose
        print("Constructing matrix...")
        flush(stdout)
    end
    L_he = MatrixOperator(
        sparse(reduce(vcat, L_row), reduce(vcat, L_col), reduce(vcat, L_val), Nado * sup_dim, Nado * sup_dim),
    )
    if verbose
        println("[DONE]")
        flush(stdout)
    end
    return M_Fermion(L_he, tier, copy(_Hsys.dimensions), Nado, sup_dim, parity, Bath, hierarchy)
end

_getBtier(M::M_Fermion) = 0
_getFtier(M::M_Fermion) = M.tier

"""
    mutable struct M_Boson <: AbstractHEOMMatrix
Heom liouvillian superoperator matrix for bosonic bath

# Fields
- `data` : the sparse matrix of HEOM liouvillian superoperator
- `tier` : the tier (cutoff) for the bath
- `dim` : the dimension of system
- `N` : the number of total ADOs
- `Nb` : the number of bosonic ADOs
- `Nf` : the number of fermionic ADOs (should be zero)
- `sup_dim` : the dimension of system superoperator
- `parity` : the parity of the density matrix (restrict to `:none` for boson)
- `bath::Vector{BosonBath}` : the vector which stores all `BosonBath` objects
- `hierarchy::HierarchyDict`: the object which contains all dictionaries for boson-bath-ADOs hierarchy.
"""
mutable struct M_Boson <: AbstractHEOMMatrix
    data::SparseMatrixCSC{ComplexF64, Int64}
    const tier::Int
    const dim::Int
    const N::Int
    const Nb::Int
    const Nf::Int
    const sup_dim::Int
    const parity::Symbol
    const bath::Vector{BosonBath}
    const hierarchy::HierarchyDict
end

function M_Boson(Hsys, tier::Int, Bath::BosonBath; verbose::Bool=true)
    return M_Boson(Hsys, tier, [Bath], verbose = verbose)
end

"""
    M_Boson(Hsys, tier, Bath; verbose=true)
Generate the boson-type Heom liouvillian superoperator matrix

# Parameters
- `Hsys` : The system Hamiltonian
- `tier::Int` : the tier (cutoff) for the bath
- `Bath::Vector{BosonBath}` : objects for different bosonic baths
- `verbose::Bool` : To display verbose output and progress bar during the process or not. Defaults to `true`.
"""
function M_Boson(        
        Hsys,
        tier::Int,
        Bath::Vector{BosonBath};
        verbose::Bool=true
    )

    # check for system dimension
    if !isValidMatrixType(Hsys)
        error("Invalid matrix \"Hsys\" (system Hamiltonian).")
    end
    Nsys,   = size(Hsys)
    sup_dim = Nsys ^ 2
    I_sup   = sparse(I, sup_dim, sup_dim)

    # the liouvillian operator for free Hamiltonian term
    Lsys = -1im * (spre(Hsys) - spost(Hsys))

    # bosonic bath
    Nado, baths, hierarchy = genBathHierarchy(Bath, tier, Nsys)
    idx2nvec = hierarchy.idx2nvec
    nvec2idx = hierarchy.nvec2idx

    # start to construct the matrix
    L_row = Int[]
    L_col = Int[]
    L_val = ComplexF64[]
    lk = SpinLock()
    if verbose
        println("Preparing block matrices for HEOM liouvillian superoperator (using $(nthreads()) threads)...")
        flush(stdout)
        prog = Progress(Nado; desc="Processing: ", PROGBAR_OPTIONS...)
    end
    @threads for idx in 1:Nado
        nvec = idx2nvec[idx]
        n_exc = sum(nvec)
        if n_exc >= 1
            sum_ω = bath_sum_ω(nvec, baths)
            op = Lsys - sum_ω * I_sup
        else
            op = Lsys
        end
        lock(lk)
        try
            add_operator!(op, L_row, L_col, L_val, Nado, idx, idx)
        finally
            unlock(lk)
        end

        count = 0
        nvec_neigh = copy(nvec)
        for bB in baths
            for k in 1:bB.Nterm
                count += 1
                n_k = nvec[count]
                if n_k >= 1
                    nvec_neigh[count] = n_k - 1
                    idx_neigh = nvec2idx[nvec_neigh]
                    
                    op = prev_grad_boson(bB, k, n_k)
                    lock(lk)
                    try
                        add_operator!(op, L_row, L_col, L_val, Nado, idx, idx_neigh)
                    finally
                        unlock(lk)
                    end
                    nvec_neigh[count] = n_k
                end
                if n_exc <= tier - 1
                    nvec_neigh[count] = n_k + 1
                    idx_neigh = nvec2idx[nvec_neigh]
                    
                    op = next_grad_boson(bB)
                    lock(lk)
                    try
                        add_operator!(op, L_row, L_col, L_val, Nado, idx, idx_neigh)
                    finally
                        unlock(lk)
                    end
                    nvec_neigh[count] = n_k
                end
            end
        end
        if verbose
            next!(prog)
        end
    end
    if verbose
        print("Constructing matrix...")
        flush(stdout)
    end
    L_he = sparse(L_row, L_col, L_val, Nado * sup_dim, Nado * sup_dim)
    if verbose
        println("[DONE]")
        flush(stdout)
    end
    return M_Boson(L_he, tier, Nsys, Nado, Nado, 0, sup_dim, :none, Bath, hierarchy)
end
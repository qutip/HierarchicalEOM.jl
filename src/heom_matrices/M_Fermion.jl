"""
    struct M_Fermion <: AbstractHEOMMatrix
Heom liouvillian superoperator matrix for fermionic bath

# Fields
- `data` : the sparse matrix of HEOM liouvillian superoperator
- `tier` : the tier (cutoff) for the bath
- `dim` : the dimension of system
- `N` : the number of total ADOs
- `sup_dim` : the dimension of system superoperator
- `parity` : the parity of the density matrix
- `bath::Vector{FermionBath}` : the vector which stores all `FermionBath` objects
- `hierarchy::HierarchyDict`: the object which contains all dictionaries for fermion-bath-ADOs hierarchy.
"""
struct M_Fermion <: AbstractHEOMMatrix
    data::SparseMatrixCSC{ComplexF64, Int64}
    tier::Int
    dim::Int
    N::Int
    sup_dim::Int
    parity::Symbol
    bath::Vector{FermionBath}
    hierarchy::HierarchyDict
end

function M_Fermion(Hsys, tier::Int, Bath::FermionBath, parity::Symbol=:even; threshold::Real=0.0, verbose::Bool=true)
    return M_Fermion(Hsys, tier, [Bath], parity, threshold= threshold, verbose = verbose)
end

"""
    M_Fermion(Hsys, tier, Bath, parity=:even; threshold=0.0, verbose=true)
Generate the fermion-type Heom liouvillian superoperator matrix

# Parameters
- `Hsys` : The system Hamiltonian
- `tier::Int` : the tier (cutoff) for the bath
- `Bath::Vector{FermionBath}` : objects for different fermionic baths
- `parity::Symbol` : The parity symbol of the density matrix (either `:odd` or `:even`). Defaults to `:even`.
- `threshold::Real` : The threshold of the importance value (see Ref. [1]). Defaults to `0.0`.
- `verbose::Bool` : To display verbose output and progress bar during the process or not. Defaults to `true`.

[1] [Phys. Rev. B 88, 235426 (2013)](https://doi.org/10.1103/PhysRevB.88.235426)
"""
function M_Fermion(        
        Hsys,
        tier::Int,
        Bath::Vector{FermionBath},
        parity::Symbol=:even;
        threshold::Real=0.0,
        verbose::Bool=true
    )

    # check parity
    if (parity != :even) && (parity != :odd)
        error("The parity symbol of density matrix should be either \":odd\" or \":even\".")
    end

    # check for system dimension
    if !isValidMatrixType(Hsys)
        error("Invalid matrix \"Hsys\" (system Hamiltonian).")
    end
    Nsys,   = size(Hsys)
    sup_dim = Nsys ^ 2
    I_sup   = sparse(I, sup_dim, sup_dim)
    
    # the liouvillian operator for free Hamiltonian term
    Lsys = -1im * (spre(Hsys) - spost(Hsys))

    # fermionic bath
    if verbose && (threshold > 0.0)
        print("Checking the importance value for each ADOs...")
        flush(stdout)
    end
    Nado, baths, hierarchy = genBathHierarchy(Bath, tier, Nsys, threshold=threshold)
    idx2nvec = hierarchy.idx2nvec
    nvec2idx = hierarchy.nvec2idx
    if verbose && (threshold > 0.0)
        println("[DONE]")
        flush(stdout)
    end

    # start to construct the matrix
    L_row = Int[]
    L_col = Int[]
    L_val = ComplexF64[]

    if verbose
        println("Preparing block matrices for HEOM liouvillian superoperator...")
        flush(stdout)
        prog = Progress(Nado; desc="Processing: ", PROGBAR_OPTIONS...)
    end
    for idx in 1:Nado
        # fermion (n tier) superoperator
        nvec = idx2nvec[idx]
        if nvec.level >= 1
            sum_γ = bath_sum_γ(nvec, baths)
            op = Lsys - sum_γ * I_sup                
        else
            op = Lsys
        end
        add_operator!(op, L_row, L_col, L_val, Nado, idx, idx)

        # fermion (n+1 & n-1 tier) superoperator
        count = 0
        nvec_neigh = copy(nvec)
        for fB in baths
            for k in 1:fB.Nterm
                count += 1
                n_k = nvec[count]

                # deal with prevous gradient
                if n_k > 0
                    Nvec_minus!(nvec_neigh, count)
                    if (threshold == 0.0) || haskey(nvec2idx, nvec_neigh)
                        idx_neigh = nvec2idx[nvec_neigh]
                        op = prev_grad_fermion(fB, k, nvec.level, sum(nvec_neigh[1:(count - 1)]), parity)
                        add_operator!(op, L_row, L_col, L_val, Nado, idx, idx_neigh)
                    end
                    Nvec_plus!(nvec_neigh, count)

                # deal with next gradient
                elseif nvec.level < tier
                    Nvec_plus!(nvec_neigh, count)
                    if (threshold == 0.0) || haskey(nvec2idx, nvec_neigh)
                        idx_neigh = nvec2idx[nvec_neigh]
                        op = next_grad_fermion(fB, nvec.level, sum(nvec_neigh[1:(count - 1)]), parity)
                        add_operator!(op, L_row, L_col, L_val, Nado, idx, idx_neigh)
                    end
                    Nvec_minus!(nvec_neigh, count)
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
    L_he = sparse(L_row, L_col, L_val, Nado * sup_dim, Nado * sup_dim)
    if verbose
        println("[DONE]")
        flush(stdout)
    end
    return M_Fermion(L_he, tier, Nsys, Nado, sup_dim, parity, Bath, hierarchy)
end
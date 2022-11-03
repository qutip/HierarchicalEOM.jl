"""
    mutable struct M_Fermion <: AbstractHEOMMatrix
Heom liouvillian superoperator matrix for fermionic bath

# Fields
- `data` : the sparse matrix of HEOM liouvillian superoperator
- `tier` : the tier (cutoff) for the bath
- `dim` : the dimension of system
- `N` : the number of total ADOs
- `Nb` : the number of bosonic ADOs (should be zero)
- `Nf` : the number of fermionic ADOs
- `sup_dim` : the dimension of system superoperator
- `parity` : the parity of the density matrix
- `bath::Vector{FermionBath}` : the vector which stores all `FermionBath` objects
- `hierarchy::HierarchyDict`: the object which contains all dictionaries for fermion-bath-ADOs hierarchy.
"""
mutable struct M_Fermion <: AbstractHEOMMatrix
    data::SparseMatrixCSC{ComplexF64, Int64}
    const tier::Int
    const dim::Int
    const N::Int
    const Nb::Int
    const Nf::Int
    const sup_dim::Int
    const parity::Symbol
    const bath::Vector{FermionBath}
    const hierarchy::HierarchyDict
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
    L_row = distribute([Int[] for _ in procs()])
    L_col = distribute([Int[] for _ in procs()])
    L_val = distribute([ComplexF64[] for _ in procs()])
    channel = RemoteChannel(() -> Channel{Bool}(), 1) # for updating the progress bar

    if verbose
        println("Preparing block matrices for HEOM liouvillian superoperator (using $(nprocs()) processors)...")
        flush(stdout)
        prog = Progress(Nado; desc="Processing: ", PROGBAR_OPTIONS...)
    end
    @sync begin # start two tasks which will be synced in the very end
        # the first task updates the progress bar
        @async while take!(channel)
            if verbose
                next!(prog)
            else
                put!(channel, false) # this tells the printing task to finish
            end
        end

        # the second task does the computation
        @async begin
            @distributed (+) for idx in 1:Nado
                nvec = idx2nvec[idx]
                n_exc = sum(nvec)
                if n_exc >= 1
                    sum_ω = bath_sum_ω(nvec, baths)
                    op = Lsys - sum_ω * I_sup                
                else
                    op = Lsys
                end
                add_operator!(op, L_row, L_col, L_val, Nado, idx, idx)

                count = 0
                nvec_neigh = copy(nvec)
                for fB in baths
                    for k in 1:fB.Nterm
                        count += 1
                        n_k = nvec[count]
                        if n_k >= 1
                            nvec_neigh[count] = n_k - 1
                            if (threshold == 0.0) || (nvec_neigh in keys(nvec2idx))
                                idx_neigh = nvec2idx[nvec_neigh]
                                op = prev_grad_fermion(fB, k, n_exc, sum(nvec_neigh[1:(count - 1)]), parity)
                                add_operator!(op, L_row, L_col, L_val, Nado, idx, idx_neigh)
                            end

                        elseif n_exc <= tier - 1
                            nvec_neigh[count] = n_k + 1
                            if (threshold == 0.0) || (nvec_neigh in keys(nvec2idx))
                                idx_neigh = nvec2idx[nvec_neigh]
                                op = next_grad_fermion(fB, n_exc, sum(nvec_neigh[1:(count - 1)]), parity)
                                add_operator!(op, L_row, L_col, L_val, Nado, idx, idx_neigh)
                            end
                            
                        else
                            continue
                        end
                        nvec_neigh[count] = n_k
                    end
                end
                if verbose
                    put!(channel, true) # trigger a progress bar update
                end
                1 # Here, returning some number 1 and reducing it somehow (+) is necessary to make the distribution happen.
            end
            put!(channel, false) # this tells the printing task to finish
        end
    end
    if verbose
        print("Constructing matrix...")
        flush(stdout)
    end
    L_he = sparse(vcat(L_row...), vcat(L_col...), vcat(L_val...), Nado * sup_dim, Nado * sup_dim)
    if verbose
        println("[DONE]")
        flush(stdout)
    end
    d_closeall()  # release all distributed arrays created from the calling process
    return M_Fermion(L_he, tier, Nsys, Nado, 0, Nado, sup_dim, parity, Bath, hierarchy)
end
"""
    mutable struct M_Boson <: AbstractHEOMMatrix
Heom liouvillian superoperator matrix for bosonic bath

# Fields
- `data` : the sparse matrix of HEOM liouvillian superoperator
- `tier` : the tier (cutoff) for the bath
- `dim` : the dimension of system
- `N` : the number of total ADOs
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
    const sup_dim::Int
    const parity::Symbol
    const bath::Vector{BosonBath}
    const hierarchy::HierarchyDict
end

function M_Boson(Hsys, tier::Int, Bath::BosonBath; threshold::Real = 0.0, verbose::Bool=true)
    return M_Boson(Hsys, tier, [Bath], threshold = threshold, verbose = verbose)
end

"""
    M_Boson(Hsys, tier, Bath; threshold=0.0, verbose=true)
Generate the boson-type Heom liouvillian superoperator matrix

# Parameters
- `Hsys` : The system Hamiltonian
- `tier::Int` : the tier (cutoff) for the bath
- `Bath::Vector{BosonBath}` : objects for different bosonic baths
- `threshold::Real` : The threshold of the importance value (see Ref. [1]). Defaults to `0.0`.
- `verbose::Bool` : To display verbose output and progress bar during the process or not. Defaults to `true`.

[1] [Phys. Rev. B 88, 235426 (2013)](https://doi.org/10.1103/PhysRevB.88.235426)
"""
function M_Boson(        
        Hsys,
        tier::Int,
        Bath::Vector{BosonBath};
        threshold::Real=0.0,
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
            try
                @distributed (+) for idx in 1:Nado
                    # boson (n tier) superoperator
                    nvec = idx2nvec[idx]
                    if nvec.level >= 1
                        sum_γ = bath_sum_γ(nvec, baths)
                        op = Lsys - sum_γ * I_sup
                    else
                        op = Lsys
                    end
                    add_operator!(op, L_row, L_col, L_val, Nado, idx, idx)

                    # boson (n+1 & n-1 tier) superoperator
                    count = 0
                    nvec_neigh = copy(nvec)
                    for bB in baths
                        for k in 1:bB.Nterm
                            count += 1
                            n_k = nvec[count]
                            
                            # deal with prevous gradient
                            if n_k > 0
                                Nvec_minus!(nvec_neigh, count)
                                if (threshold == 0.0) || haskey(nvec2idx, nvec_neigh)
                                    idx_neigh = nvec2idx[nvec_neigh]
                                    op = prev_grad_boson(bB, k, n_k)
                                    add_operator!(op, L_row, L_col, L_val, Nado, idx, idx_neigh)
                                end
                                Nvec_plus!(nvec_neigh, count)
                            end

                            # deal with next gradient
                            if nvec.level < tier
                                Nvec_plus!(nvec_neigh, count)
                                if (threshold == 0.0) || haskey(nvec2idx, nvec_neigh)
                                    idx_neigh = nvec2idx[nvec_neigh]
                                    op = next_grad_boson(bB)
                                    add_operator!(op, L_row, L_col, L_val, Nado, idx, idx_neigh)
                                end
                                Nvec_minus!(nvec_neigh, count)
                            end
                        end
                    end
                    if verbose
                        put!(channel, true) # trigger a progress bar update
                    end
                    1 # Here, returning some number 1 and reducing it somehow (+) is necessary to make the distribution happen.
                end
                put!(channel, false) # this tells the printing task to finish

            catch e
                put!(channel, false) # this tells the printing task to finish
                throw(e)
            end
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
    return M_Boson(L_he, tier, Nsys, Nado, sup_dim, :none, Bath, hierarchy)
end
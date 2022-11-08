"""
    mutable struct M_Boson_Fermion <: AbstractHEOMMatrix
Heom liouvillian superoperator matrix for mixtured (bosonic and fermionic) bath 

# Fields
- `data` : the sparse matrix of HEOM liouvillian superoperator
- `Btier` : the tier (cutoff) for bosonic bath
- `Ftier` : the tier (cutoff) for fermionic bath
- `dim` : the dimension of system
- `N` : the number of total ADOs
- `sup_dim` : the dimension of system superoperator
- `parity` : the parity of the density matrix
- `Bbath::Vector{BosonBath}` : the vector which stores all `BosonBath` objects
- `Fbath::Vector{FermionBath}` : the vector which stores all `FermionBath` objects
- `hierarchy::MixHierarchyDict`: the object which contains all dictionaries for mixed-bath-ADOs hierarchy.
"""
mutable struct M_Boson_Fermion <: AbstractHEOMMatrix
    data::SparseMatrixCSC{ComplexF64, Int64}
    const Btier::Int
    const Ftier::Int
    const dim::Int
    const N::Int
    const sup_dim::Int
    const parity::Symbol
    const Bbath::Vector{BosonBath}
    const Fbath::Vector{FermionBath}
    const hierarchy::MixHierarchyDict
end

function M_Boson_Fermion(Hsys, Btier::Int, Ftier::Int, Bbath::BosonBath, Fbath::FermionBath, parity::Symbol=:even; threshold::Real = 0.0, verbose::Bool=true)
    return M_Boson_Fermion(Hsys, Btier, Ftier, [Bbath], [Fbath], parity, threshold = threshold, verbose = verbose)
end

function M_Boson_Fermion(Hsys, Btier::Int, Ftier::Int, Bbath::Vector{BosonBath}, Fbath::FermionBath, parity::Symbol=:even; threshold::Real = 0.0, verbose::Bool=true)
    return M_Boson_Fermion(Hsys, Btier, Ftier, Bbath, [Fbath], parity, threshold = threshold, verbose = verbose)
end

function M_Boson_Fermion(Hsys, Btier::Int, Ftier::Int, Bbath::BosonBath, Fbath::Vector{FermionBath}, parity::Symbol=:even; threshold::Real = 0.0, verbose::Bool=true)
    return M_Boson_Fermion(Hsys, Btier, Ftier, [Bbath], Fbath, parity, threshold = threshold, verbose = verbose)
end

"""
    M_Boson_Fermion(Hsys, Btier, Ftier, Bbath, Fbath, parity=:even; threshold=0.0, verbose=true)
Generate the boson-fermion-type Heom liouvillian superoperator matrix

# Parameters
- `Hsys` : The system Hamiltonian
- `Btier::Int` : the tier (cutoff) for the bosonic bath
- `Ftier::Int` : the tier (cutoff) for the fermionic bath
- `Bbath::Vector{BosonBath}` : objects for different bosonic baths
- `Fbath::Vector{FermionBath}` : objects for different fermionic baths
- `parity::Symbol` : The parity symbol of the density matrix (either `:odd` or `:even`). Defaults to `:even`.
- `threshold::Real` : The threshold of the importance value (see Ref. [1, 2]). Defaults to `0.0`.
- `verbose::Bool` : To display verbose output and progress bar during the process or not. Defaults to `true`.

[1] [Phys. Rev. B  88, 235426 (2013)](https://doi.org/10.1103/PhysRevB.88.235426)
[2] [Phys. Rev. B 103, 235413 (2021)](https://doi.org/10.1103/PhysRevB.103.235413)
"""
function M_Boson_Fermion(        
        Hsys,
        Btier::Int,
        Ftier::Int,
        Bbath::Vector{BosonBath},
        Fbath::Vector{FermionBath},
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

    # check for bosonic and fermionic bath
    if verbose && (threshold > 0.0)
        print("Checking the importance value for each ADOs...")
        flush(stdout)
    end
    Nado, baths_b, baths_f, hierarchy = genBathHierarchy(Bbath, Fbath, Btier, Ftier, Nsys, threshold = threshold)
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
                # boson and fermion (n tier) superoperator
                sum_γ   = 0.0
                nvec_b, nvec_f = idx2nvec[idx]
                if nvec_b.level >= 1
                    sum_γ += bath_sum_γ(nvec_b, baths_b)
                end
                if nvec_f.level >= 1
                    sum_γ += bath_sum_γ(nvec_f, baths_f)
                end
                add_operator!(Lsys - sum_γ * I_sup, L_row, L_col, L_val, Nado, idx, idx)
                
                # boson (n+1 & n-1 tier) superoperator
                count = 0
                nvec_neigh = copy(nvec_b)
                for bB in baths_b
                    for k in 1:bB.Nterm
                        count += 1
                        n_k = nvec_b[count]

                        # deal with prevous gradient
                        if n_k > 0
                            Nvec_minus!(nvec_neigh, count)
                            if (threshold == 0.0) || haskey(nvec2idx, (nvec_neigh, nvec_f))
                                idx_neigh = nvec2idx[(nvec_neigh, nvec_f)]
                                op = prev_grad_boson(bB, k, n_k)
                                add_operator!(op, L_row, L_col, L_val, Nado, idx, idx_neigh)
                            end
                            Nvec_plus!(nvec_neigh, count)
                        end

                        # deal with next gradient
                        if nvec_b.level < Btier
                            Nvec_plus!(nvec_neigh, count)
                            if (threshold == 0.0) || haskey(nvec2idx, (nvec_neigh, nvec_f))
                                idx_neigh = nvec2idx[(nvec_neigh, nvec_f)]
                                op = next_grad_boson(bB)
                                add_operator!(op, L_row, L_col, L_val, Nado, idx, idx_neigh)
                            end
                            Nvec_minus!(nvec_neigh, count)
                        end
                    end
                end
                
                # fermion (n+1 & n-1 tier) superoperator
                count = 0
                nvec_neigh = copy(nvec_f)
                for fB in baths_f
                    for k in 1:fB.Nterm
                        count += 1
                        n_k = nvec_f[count]

                        # deal with prevous gradient
                        if n_k > 0
                            Nvec_minus!(nvec_neigh, count)
                            if (threshold == 0.0) || haskey(nvec2idx, (nvec_b, nvec_neigh))
                                idx_neigh = nvec2idx[(nvec_b, nvec_neigh)]
                                op = prev_grad_fermion(fB, k, nvec_f.level, sum(nvec_neigh[1:(count - 1)]), parity)
                                add_operator!(op, L_row, L_col, L_val, Nado, idx, idx_neigh)
                            end
                            Nvec_plus!(nvec_neigh, count)

                        # deal with next gradient
                        elseif nvec_f.level < Ftier
                            Nvec_plus!(nvec_neigh, count)
                            if (threshold == 0.0) || haskey(nvec2idx, (nvec_b, nvec_neigh))
                                idx_neigh = nvec2idx[(nvec_b, nvec_neigh)]
                                op = next_grad_fermion(fB, nvec_f.level, sum(nvec_neigh[1:(count - 1)]), parity)
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
    return M_Boson_Fermion(L_he, Btier, Ftier, Nsys, Nado, sup_dim, parity, Bbath, Fbath, hierarchy)
end
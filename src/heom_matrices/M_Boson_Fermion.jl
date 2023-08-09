@doc raw"""
    struct M_Boson_Fermion <: AbstractHEOMLSMatrix
HEOM Liouvillian superoperator matrix for mixtured (bosonic and fermionic) bath 

# Fields
- `data` : the sparse matrix of HEOM Liouvillian superoperator
- `Btier` : the tier (cutoff level) for bosonic hierarchy
- `Ftier` : the tier (cutoff level) for fermionic hierarchy
- `dim` : the dimension of system
- `N` : the number of total ADOs
- `sup_dim` : the dimension of system superoperator
- `parity` : the parity label of the fermionic system (usually `:even`, only set as `:odd` for calculating spectrum of fermionic system).
- `Bbath::Vector{BosonBath}` : the vector which stores all `BosonBath` objects
- `Fbath::Vector{FermionBath}` : the vector which stores all `FermionBath` objects
- `hierarchy::MixHierarchyDict`: the object which contains all dictionaries for mixed-bath-ADOs hierarchy.
"""
struct M_Boson_Fermion <: AbstractHEOMLSMatrix
    data::SparseMatrixCSC{ComplexF64, Int64}
    Btier::Int
    Ftier::Int
    dim::Int
    N::Int
    sup_dim::Int
    parity::Symbol
    Bbath::Vector{BosonBath}
    Fbath::Vector{FermionBath}
    hierarchy::MixHierarchyDict
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

@doc raw"""
    M_Boson_Fermion(Hsys, Btier, Ftier, Bbath, Fbath, parity=:even; threshold=0.0, verbose=true)
Generate the boson-fermion-type HEOM Liouvillian superoperator matrix

# Parameters
- `Hsys` : The time-independent system Hamiltonian
- `Btier::Int` : the tier (cutoff level) for the bosonic bath
- `Ftier::Int` : the tier (cutoff level) for the fermionic bath
- `Bbath::Vector{BosonBath}` : objects for different bosonic baths
- `Fbath::Vector{FermionBath}` : objects for different fermionic baths
- `parity::Symbol` : the parity label of the fermionic system (only set as `:odd` for calculating spectrum of fermionic system). Defaults to `:even`.
- `threshold::Real` : The threshold of the importance value (see Ref. [1, 2]). Defaults to `0.0`.
- `verbose::Bool` : To display verbose output and progress bar during the process or not. Defaults to `true`.

Note that the parity only need to be set as `:odd` when the system contains fermion systems and you need to calculate the spectrum of it.

[1] [Phys. Rev. B  88, 235426 (2013)](https://doi.org/10.1103/PhysRevB.88.235426)
[2] [Phys. Rev. B 103, 235413 (2021)](https://doi.org/10.1103/PhysRevB.103.235413)
"""
@noinline function M_Boson_Fermion(        
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
        error("The parity symbol of density matrix should be either \":even\" or \":odd\".")
    end

    # check for system dimension
    if !isValidMatrixType(Hsys)
        error("Invalid matrix \"Hsys\" (system Hamiltonian).")
    end
    Nsys,   = size(Hsys)
    sup_dim = Nsys ^ 2
    I_sup   = sparse(I, sup_dim, sup_dim)

    # the Liouvillian operator for free Hamiltonian term
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
    Nthread = nthreads()
    L_row = [Int[] for _ in 1:Nthread]
    L_col = [Int[] for _ in 1:Nthread]
    L_val = [ComplexF64[] for _ in 1:Nthread]

    if verbose
        println("Preparing block matrices for HEOM Liouvillian superoperator (using $(Nthread) threads)...")
        flush(stdout)
        prog = Progress(Nado; desc="Processing: ", PROGBAR_OPTIONS...)
    end
    @threads for idx in 1:Nado
        tID = threadid()

        # boson and fermion (current level) superoperator
        sum_γ   = 0.0
        nvec_b, nvec_f = idx2nvec[idx]
        if nvec_b.level >= 1
            sum_γ += bath_sum_γ(nvec_b, baths_b)
        end
        if nvec_f.level >= 1
            sum_γ += bath_sum_γ(nvec_f, baths_f)
        end
        add_operator!(Lsys - sum_γ * I_sup, L_row[tID], L_col[tID], L_val[tID], Nado, idx, idx)
        
        # connect to bosonic (n+1)th- & (n-1)th- level superoperator
        mode = 0
        nvec_neigh = copy(nvec_b)
        for bB in baths_b
            for k in 1:bB.Nterm
                mode += 1
                n_k = nvec_b[mode]

                # connect to bosonic (n-1)th-level superoperator
                if n_k > 0
                    Nvec_minus!(nvec_neigh, mode)
                    if (threshold == 0.0) || haskey(nvec2idx, (nvec_neigh, nvec_f))
                        idx_neigh = nvec2idx[(nvec_neigh, nvec_f)]
                        op = _D_op(bB, k, n_k)
                        add_operator!(op, L_row[tID], L_col[tID], L_val[tID], Nado, idx, idx_neigh)
                    end
                    Nvec_plus!(nvec_neigh, mode)
                end

                # connect to bosonic (n+1)th-level superoperator
                if nvec_b.level < Btier
                    Nvec_plus!(nvec_neigh, mode)
                    if (threshold == 0.0) || haskey(nvec2idx, (nvec_neigh, nvec_f))
                        idx_neigh = nvec2idx[(nvec_neigh, nvec_f)]
                        op = _B_op(bB)
                        add_operator!(op, L_row[tID], L_col[tID], L_val[tID], Nado, idx, idx_neigh)
                    end
                    Nvec_minus!(nvec_neigh, mode)
                end
            end
        end
        
        # connect to fermionic (n+1)th- & (n-1)th- level superoperator
        mode = 0
        nvec_neigh = copy(nvec_f)
        for fB in baths_f
            for k in 1:fB.Nterm
                mode += 1
                n_k = nvec_f[mode]

                # connect to fermionic (n-1)th-level superoperator
                if n_k > 0
                    Nvec_minus!(nvec_neigh, mode)
                    if (threshold == 0.0) || haskey(nvec2idx, (nvec_b, nvec_neigh))
                        idx_neigh = nvec2idx[(nvec_b, nvec_neigh)]
                        op = _C_op(fB, k, nvec_f.level, sum(nvec_neigh[1:(mode - 1)]), parity)
                        add_operator!(op, L_row[tID], L_col[tID], L_val[tID], Nado, idx, idx_neigh)
                    end
                    Nvec_plus!(nvec_neigh, mode)

                # connect to fermionic (n+1)th-level superoperator
                elseif nvec_f.level < Ftier
                    Nvec_plus!(nvec_neigh, mode)
                    if (threshold == 0.0) || haskey(nvec2idx, (nvec_b, nvec_neigh))
                        idx_neigh = nvec2idx[(nvec_b, nvec_neigh)]
                        op = _A_op(fB, nvec_f.level, sum(nvec_neigh[1:(mode - 1)]), parity)
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
    L_he = sparse(reduce(vcat, L_row), reduce(vcat, L_col), reduce(vcat, L_val), Nado * sup_dim, Nado * sup_dim)
    if verbose 
        println("[DONE]") 
        flush(stdout)
    end
    return M_Boson_Fermion(L_he, Btier, Ftier, Nsys, Nado, sup_dim, parity, Bbath, Fbath, hierarchy)
end
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
- `baths::Vector{FermionBath}` : the vector which stores all `FermionBath` objects
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
    const baths::Vector{FermionBath}
    const hierarchy::HierarchyDict
end

function M_Fermion(Hsys, tier::Int, Bath::FermionBath, parity::Symbol=:even; progressBar::Bool=true)
    return M_Fermion(Hsys, tier, [Bath], parity, progressBar = progressBar)
end

"""
    M_Fermion(Hsys, tier, Bath, parity=:even; progressBar=true)
Generate the fermion-type Heom liouvillian superoperator matrix

# Parameters
- `Hsys` : The system Hamiltonian
- `tier::Int` : the tier (cutoff) for the bath
- `Bath::Vector{FermionBath}` : objects for different fermionic baths
- `parity::Symbol` : The parity symbol of the density matrix (either `:odd` or `:even`). Defaults to `:even`.
- `progressBar::Bool` : Display progress bar during the process or not. Defaults to `true`.
"""
function M_Fermion(        
        Hsys,
        tier::Int,
        Bath::Vector{FermionBath},
        parity::Symbol=:even;
        progressBar::Bool=true
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
    Nado, bath, hierarchy = genBathHierarchy(Bath, tier, Nsys)
    idx2ado = hierarchy.idx2ado
    ado2idx = hierarchy.ado2idx

    # start to construct the matrix
    L_row = distribute([Int[] for _ in procs()])
    L_col = distribute([Int[] for _ in procs()])
    L_val = distribute([ComplexF64[] for _ in procs()])
    channel = RemoteChannel(() -> Channel{Bool}(), 1) # for updating the progress bar

    println("Preparing block matrices for HEOM liouvillian superoperator (using $(nprocs()) processors)...")
    if progressBar
        prog = Progress(Nado; desc="Processing: ", PROGBAR_OPTIONS...)
    else
        println("Processing...")
        flush(stdout)
    end
    @sync begin # start two tasks which will be synced in the very end
        # the first task updates the progress bar
        @async while take!(channel)
            if progressBar
                next!(prog)
            else
                put!(channel, false) # this tells the printing task to finish
            end
        end

        # the second task does the computation
        @async begin
            @distributed (+) for idx in 1:Nado
                ado = idx2ado[idx]
                n_exc = sum(ado)
                if n_exc >= 1
                    sum_ω = bath_sum_ω(ado, baths)
                    op = Lsys - sum_ω * I_sup                
                else
                    op = Lsys
                end
                add_operator!(op, L_row, L_col, L_val, Nado, idx, idx)

                count = 0
                ado_neigh = copy(ado)
                for fB in bath
                    for k in 1:fB.Nterm
                        count += 1
                        n_k = ado[count]
                        if n_k >= 1
                            ado_neigh[count] = n_k - 1
                            idx_neigh = ado2idx[ado_neigh]
                            op = prev_grad_fermion(fB, k, n_exc, sum(ado_neigh[1:(count - 1)]), parity)

                        elseif n_exc <= tier - 1
                            ado_neigh[count] = n_k + 1
                            idx_neigh = ado2idx[ado_neigh]
                            op = next_grad_fermion(fB, n_exc, sum(ado_neigh[1:(count - 1)]), parity)
                        
                        else
                            continue
                        end
                        add_operator!(op, L_row, L_col, L_val, Nado, idx, idx_neigh)
                        
                        ado_neigh[count] = n_k
                    end
                end
                if progressBar
                    put!(channel, true) # trigger a progress bar update
                end
                1 # Here, returning some number 1 and reducing it somehow (+) is necessary to make the distribution happen.
            end
            put!(channel, false) # this tells the printing task to finish
        end
    end
    print("Constructing matrix...")
    flush(stdout)
    L_he = sparse(vcat(L_row...), vcat(L_col...), vcat(L_val...), Nado * sup_dim, Nado * sup_dim)
    println("[DONE]")

    return M_Fermion(L_he, tier, Nsys, Nado, 0, Nado, sup_dim, parity, Bath, hierarchy)
end
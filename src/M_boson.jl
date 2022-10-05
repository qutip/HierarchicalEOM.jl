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
    idx2ado = hierarchy.idx2ado
    ado2idx = hierarchy.ado2idx

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
                for bB in baths
                    for k in 1:bB.Nterm
                        count += 1
                        n_k = ado[count]
                        if n_k >= 1
                            ado_neigh[count] = n_k - 1
                            idx_neigh = ado2idx[ado_neigh]
                            
                            op = prev_grad_boson(bB, k, n_k)
                            add_operator!(op, L_row, L_col, L_val, Nado, idx, idx_neigh)

                            ado_neigh[count] = n_k
                        end
                        if n_exc <= tier - 1
                            ado_neigh[count] = n_k + 1
                            idx_neigh = ado2idx[ado_neigh]
                            
                            op = next_grad_boson(bB)
                            add_operator!(op, L_row, L_col, L_val, Nado, idx, idx_neigh)
                            
                            ado_neigh[count] = n_k
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
    
    return M_Boson(L_he, tier, Nsys, Nado, Nado, 0, sup_dim, :none, Bath, hierarchy)
end
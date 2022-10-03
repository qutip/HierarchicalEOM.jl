"""
    struct M_Fermion <: AbstractHEOMMatrix
Heom liouvillian superoperator matrix for fermionic bath

# Fields
- `data` : the sparse matrix of HEOM liouvillian superoperator
- `tier` : the tier (cutoff) for the bath
- `dim` : the dimension of system
- `N` : the number of total states
- `Nb` : the number of bosonic states (should be zero)
- `Nf` : the number of fermionic states
- `sup_dim` : the dimension of system superoperator
- `parity` : the parity of the density matrix
- `ado2idx` : the ADO-to-index dictionary
"""
struct M_Fermion <: AbstractHEOMMatrix
    data::SparseMatrixCSC{ComplexF64, Int64}
    tier::Int
    dim::Int
    N::Int
    Nb::Int
    Nf::Int
    sup_dim::Int
    parity::Symbol
    ado2idx::OrderedDict{Vector{Int}, Int}
end

function M_Fermion(Hsys, tier::Int, Bath::FermionBath, parity::Symbol=:even; progressBar::Bool=true)
    return M_Fermion(Hsys, tier, [Bath], parity, progressBar = progressBar)
end

"""
    M_Fermion(Hsys, tier, Bath, parity=:even; progressBar=true)
Generate the fermion-type Heom matrix

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

    # check for fermionic bath
    if length(Bath) > 1
        baths = CombinedBath(Nsys, Bath)
    else
        baths = Bath[1]
    end
    bath       = baths.bath
    N_exp_term = baths.Nterm

    # get ADOs dictionary
    N_he, ado2idx_ordered, idx2ado = ADOs_dictionary(fill(2, N_exp_term), tier)
    ado2idx = Dict(ado2idx_ordered)

    # start to construct the matrix
    L_row = distribute([Int[] for _ in procs()])
    L_col = distribute([Int[] for _ in procs()])
    L_val = distribute([ComplexF64[] for _ in procs()])
    channel = RemoteChannel(() -> Channel{Bool}(), 1) # for updating the progress bar

    println("Preparing block matrices for HEOM liouvillian superoperator (using $(nprocs()) processors)...")
    if progressBar
        prog = Progress(N_he; desc="Processing: ", PROGBAR_OPTIONS...)
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
            @distributed (+) for idx in 1:N_he
                state = idx2ado[idx]
                n_exc = sum(state)
                if n_exc >= 1
                    sum_ω = bath_sum_ω(state, baths)
                    op = Lsys - sum_ω * I_sup                
                else
                    op = Lsys
                end
                add_operator!(op, L_row, L_col, L_val, N_he, idx, idx)

                count = 0
                state_neigh = copy(state)
                for fB in bath
                    for k in 1:fB.Nterm
                        count += 1
                        n_k = state[count]
                        if n_k >= 1
                            state_neigh[count] = n_k - 1
                            idx_neigh = ado2idx[state_neigh]
                            op = prev_grad_fermion(fB, k, n_exc, sum(state_neigh[1:(count - 1)]), parity)

                        elseif n_exc <= tier - 1
                            state_neigh[count] = n_k + 1
                            idx_neigh = ado2idx[state_neigh]
                            op = next_grad_fermion(fB, n_exc, sum(state_neigh[1:(count - 1)]), parity)
                        
                        else
                            continue
                        end
                        add_operator!(op, L_row, L_col, L_val, N_he, idx, idx_neigh)
                        
                        state_neigh[count] = n_k
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
    L_he = sparse(vcat(L_row...), vcat(L_col...), vcat(L_val...), N_he * sup_dim, N_he * sup_dim)
    println("[DONE]")

    return M_Fermion(L_he, tier, Nsys, N_he, 0, N_he, sup_dim, parity, ado2idx_ordered)
end
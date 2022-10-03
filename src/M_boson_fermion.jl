"""
    mutable struct M_Boson_Fermion <: AbstractHEOMMatrix
Heom liouvillian superoperator matrix for mixtured (bosonic and fermionic) bath 

# Fields
- `data` : the sparse matrix of HEOM liouvillian superoperator
- `tier_b` : the tier (cutoff) for bosonic bath
- `tier_f` : the tier (cutoff) for fermionic bath
- `dim` : the dimension of system
- `N` : the number of total ADOs
- `Nb` : the number of bosonic ADOs
- `Nf` : the number of fermionic ADOs
- `sup_dim` : the dimension of system superoperator
- `parity` : the parity of the density matrix
- `ado2idx_b` : the bosonic ADO-to-index dictionary
- `ado2idx_f` : the fermionic ADO-to-index dictionary
"""
mutable struct M_Boson_Fermion <: AbstractHEOMMatrix
    data::SparseMatrixCSC{ComplexF64, Int64}
    const tier_b::Int
    const tier_f::Int
    const dim::Int
    const N::Int
    const Nb::Int
    const Nf::Int
    const sup_dim::Int
    const parity::Symbol
    const ado2idx_b::OrderedDict{Vector{Int}, Int}
    const ado2idx_f::OrderedDict{Vector{Int}, Int}
end

function M_Boson_Fermion(Hsys, tier_b::Int, tier_f::Int, Bath_b::BosonBath, Bath_f::FermionBath, parity::Symbol=:even; progressBar::Bool=true)
    return M_Boson_Fermion(Hsys, tier_b, tier_f, [Bath_b], [Bath_f], parity, progressBar = progressBar)
end

function M_Boson_Fermion(Hsys, tier_b::Int, tier_f::Int, Bath_b::Vector{BosonBath}, Bath_f::FermionBath, parity::Symbol=:even; progressBar::Bool=true)
    return M_Boson_Fermion(Hsys, tier_b, tier_f, Bath_b, [Bath_f], parity, progressBar = progressBar)
end

function M_Boson_Fermion(Hsys, tier_b::Int, tier_f::Int, Bath_b::BosonBath, Bath_f::Vector{FermionBath}, parity::Symbol=:even; progressBar::Bool=true)
    return M_Boson_Fermion(Hsys, tier_b, tier_f, [Bath_b], Bath_f, parity, progressBar = progressBar)
end

"""
    M_Boson_Fermion(Hsys, tier_b, tier_f, Bath_b, Bath_f, parity=:even; progressBar=true)
Generate the boson-fermion-type Heom matrix

# Parameters
- `Hsys` : The system Hamiltonian
- `tier_b::Int` : the tier (cutoff) for the bosonic bath
- `tier_f::Int` : the tier (cutoff) for the fermionic bath
- `Bath_b::Vector{BosonBath}` : objects for different bosonic baths
- `Bath_f::Vector{FermionBath}` : objects for different fermionic baths
- `parity::Symbol` : The parity symbol of the density matrix (either `:odd` or `:even`). Defaults to `:even`.
- `progressBar::Bool` : Display progress bar during the process or not. Defaults to `true`.
"""
function M_Boson_Fermion(        
        Hsys,
        tier_b::Int,
        tier_f::Int,
        Bath_b::Vector{BosonBath},
        Bath_f::Vector{FermionBath},
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

    # check for bosonic bath
    if length(Bath_b) > 1
        baths_b = CombinedBath(Nsys, Bath_b)
    else
        baths_b = Bath_b[1]
    end
    bath_b       = baths_b.bath
    N_exp_term_b = baths_b.Nterm

    # check for fermionic bath
    if length(Bath_f) > 1
        baths_f = CombinedBath(Nsys, Bath_f)
    else
        baths_f = Bath_f[1]
    end
    bath_f       = baths_f.bath
    N_exp_term_f = baths_f.Nterm

    # get ADOs dictionary
    Nado_b, ado2idx_b_ordered, idx2ado_b = ADOs_dictionary(fill((tier_b + 1), N_exp_term_b), tier_b)
    Nado_f, ado2idx_f_ordered, idx2ado_f = ADOs_dictionary(fill(2, N_exp_term_f), tier_f)
    Nado_tot = Nado_b * Nado_f
    ado2idx_b = Dict(ado2idx_b_ordered)
    ado2idx_f = Dict(ado2idx_f_ordered)

    # start to construct the matrix
    L_row = distribute([Int[] for _ in procs()])
    L_col = distribute([Int[] for _ in procs()])
    L_val = distribute([ComplexF64[] for _ in procs()])
    channel = RemoteChannel(() -> Channel{Bool}(), 1) # for updating the progress bar

    println("Preparing block matrices for HEOM liouvillian superoperator (using $(nprocs()) processors)...")
    if progressBar
        prog = Progress(Nado_b + Nado_f; desc="Processing: ", PROGBAR_OPTIONS...)
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
            @distributed (+) for idx_b in 1:Nado_b
                # diagonal (boson)
                sum_ω   = 0.0
                ado_b = idx2ado_b[idx_b]
                n_exc_b = sum(ado_b) 
                idx = (idx_b - 1) * Nado_f
                if n_exc_b >= 1
                    sum_ω += bath_sum_ω(ado_b, baths_b)
                end

                # diagonal (fermion)
                for idx_f in 1:Nado_f
                    ado_f = idx2ado_f[idx_f]
                    n_exc_f = sum(ado_f)
                    if n_exc_f >= 1
                        sum_ω += bath_sum_ω(ado_f, baths_f)
                    end
                    add_operator!(Lsys - sum_ω * I_sup, L_row, L_col, L_val, Nado_tot, idx + idx_f, idx + idx_f)
                end
                
                # off-diagonal (boson)
                count = 0
                ado_neigh = copy(ado_b)
                for bB in bath_b
                    for k in 1:bB.Nterm
                        count += 1
                        n_k = ado_b[count]
                        if n_k >= 1
                            ado_neigh[count] = n_k - 1
                            idx_neigh = ado2idx_b[ado_neigh]
                            
                            op = prev_grad_boson(bB, k, n_k)
                            for idx_f in 1:Nado_f
                                add_operator!(op, L_row, L_col, L_val, Nado_tot, (idx + idx_f), (idx_neigh - 1) * Nado_f + idx_f)
                            end
                            
                            ado_neigh[count] = n_k
                        end
                        if n_exc_b <= tier_b - 1
                            ado_neigh[count] = n_k + 1
                            idx_neigh = ado2idx_b[ado_neigh]
                            
                            op = next_grad_boson(bB)
                            for idx_f in 1:Nado_f
                                add_operator!(op, L_row, L_col, L_val, Nado_tot, (idx + idx_f), (idx_neigh - 1) * Nado_f + idx_f)
                            end

                            ado_neigh[count] = n_k
                        end
                    end
                end
                if progressBar
                    put!(channel, true) # trigger a progress bar update
                end
                1 # Here, returning some number 1 and reducing it somehow (+) is necessary to make the distribution happen.
            end

            # fermion (n+1 & n-1 tier) superoperator
            @distributed (+) for idx_f in 1:Nado_f
                ado_f = idx2ado_f[idx_f]
                n_exc_f = sum(ado_f)

                count = 0
                ado_neigh = copy(ado_f)
                for fB in bath_f
                    for k in 1:fB.Nterm
                        count += 1
                        n_k = ado_f[count]
                        if n_k >= 1
                            ado_neigh[count] = n_k - 1
                            idx_neigh = ado2idx_f[ado_neigh]
                            op = prev_grad_fermion(fB, k, n_exc_f, sum(ado_neigh[1:(count - 1)]), parity)

                        elseif n_exc_f <= tier_f - 1
                            ado_neigh[count] = n_k + 1
                            idx_neigh = ado2idx_f[ado_neigh]
                            op = next_grad_fermion(fB, n_exc_f, sum(ado_neigh[1:(count - 1)]), parity)

                        else
                            continue
                        end

                        for idx_b in 1:Nado_b
                            idx = (idx_b - 1) * Nado_f
                            add_operator!(op, L_row, L_col, L_val, Nado_tot, idx + idx_f, idx + idx_neigh)
                        end

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
    L_he = sparse(vcat(L_row...), vcat(L_col...), vcat(L_val...), Nado_tot * sup_dim, Nado_tot * sup_dim)
    println("[DONE]")

    return M_Boson_Fermion(L_he, tier_b, tier_f, Nsys, Nado_tot, Nado_b, Nado_f, sup_dim, parity, ado2idx_b_ordered, ado2idx_f_ordered)
end
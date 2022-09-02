"""
# `M_Boson_Fermion <: AbstractHEOMMatrix`
Heom matrix for mixtured (bosonic and fermionic) bath 

## Fields
- `data::SparseMatrixCSC{ComplexF64, Int64}` : the sparse matrix
- `tier_b::Int`   : the tier (cutoff) for bosonic bath
- `tier_f::Int`   : the tier (cutoff) for fermionic bath
- `dim::Int`  : the dimension of system
- `N::Int`   : the number of total states
- `Nb::Int`   : the number of bosonic states
- `Nf::Int`   : the number of fermionic states
- `sup_dim::Int`: the dimension of system superoperator
- `parity::Symbol`: the parity of the density matrix
- `ado2idx_b::OrderedDict{Vector{Int}, Int}`: the bosonic ADO-to-index dictionary
- `ado2idx_f::OrderedDict{Vector{Int}, Int}`: the fermionic ADO-to-index dictionary

## Constructor
`M_Boson_Fermion(Hsys, tier_b, tier_f, bath_b, bath_f, parity; [progressBar])`

- `Hsys::AbstractMatrix` : The system Hamiltonian
- `tier_b::Int` : the tier (cutoff) for the bosonic bath
- `tier_f::Int` : the tier (cutoff) for the fermionic bath
- `bath_b::BosonicBath` : an object for the bosonic bath correlation
- `bath_f::FermionicBath` : an object for the fermionic bath correlation
- `parity::Symbol` : The parity symbol of the density matrix (either `:odd` or `:even`). Defaults to `:even`.
- `progressBar::Bool` : Display progress bar during the process or not. Defaults to `true`.
"""
mutable struct M_Boson_Fermion <: AbstractHEOMMatrix
    data::SparseMatrixCSC{ComplexF64, Int64}
    tier_b::Int
    tier_f::Int
    dim::Int
    N::Int
    Nb::Int
    Nf::Int
    sup_dim::Int
    parity::Symbol
    ado2idx_b::OrderedDict{Vector{Int}, Int}
    ado2idx_f::OrderedDict{Vector{Int}, Int}
    
    function M_Boson_Fermion(        
            Hsys::AbstractMatrix,
            tier_b::Int,
            tier_f::Int,
            bath_b::BosonicBath,
            bath_f::FermionicBath,
            parity::Symbol=:even;
            progressBar::Bool=true
        )

        if (parity != :even) || (parity != :odd)
            error("The parity symbol of density matrix should be either \":odd\" or \":even\".")
        end

        Nsys,   = size(Hsys)
        sup_dim = Nsys ^ 2
        I_sup   = sparse(I, sup_dim, sup_dim)

        c_list = bath_b.c_list
        ν_list = bath_b.ν_list
        η_list = bath_f.η_list
        γ_list = bath_f.γ_list
        Coup_Op_b = bath_b.coupOP
        Coup_Op_f = bath_f.coupOP
        N_oper_f  = bath_f.N_oper
        N_exp_term_b = bath_b.N_term
        N_exp_term_f = bath_f.N_term

        dims_b    = [(tier_b + 1) for i in 1:N_exp_term_b]
        dims_f    = [2 for i in 1:(N_exp_term_f * N_oper_f)]
    
        spreQ_b   = spre(Coup_Op_b)
        spostQ_b  = spost(Coup_Op_b)
        commQ_b   = spreQ_b - spostQ_b
        spreQ_f   = spre.(Coup_Op_f)
        spostQ_f  = spost.(Coup_Op_f)
        spreQd_f  = spre.(adjoint.(Coup_Op_f))
        spostQd_f = spost.(adjoint.(Coup_Op_f))

        # get ADOs dictionary
        N_he_b, ado2idx_b_ordered, idx2ado_b = ADOs_dictionary(dims_b, tier_b)
        N_he_f, ado2idx_f_ordered, idx2ado_f = ADOs_dictionary(dims_f, tier_f)
        N_he_tot = N_he_b * N_he_f
        ado2idx_b = Dict(ado2idx_b_ordered)
        ado2idx_f = Dict(ado2idx_f_ordered)

        # start to construct the matrix
        L_row = distribute([Int[] for _ in procs()])
        L_col = distribute([Int[] for _ in procs()])
        L_val = distribute([ComplexF64[] for _ in procs()])
        channel = RemoteChannel(() -> Channel{Bool}(), 1) # for updating the progress bar

        println("Start constructing hierarchy matrix...(using $(nprocs()) processors)")
        if progressBar
            prog = Progress(N_he_b + N_he_f; desc="Processing: ", PROGBAR_OPTIONS...)
        else
            println("Processing...")
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
                @distributed (+) for idx_b in 1:N_he_b
                    # diagonal (boson)
                    sum_ω   = 0.0
                    state_b = idx2ado_b[idx_b]
                    n_exc_b = sum(state_b) 
                    idx = (idx_b - 1) * N_he_f
                    if n_exc_b >= 1
                        for k in 1:N_exp_term_b
                            if state_b[k] > 0
                                sum_ω += state_b[k] * ν_list[k]
                            end
                        end
                    end

                    # diagonal (fermion)
                    for idx_f in 1:N_he_f
                        state_f = idx2ado_f[idx_f]
                        n_exc_f = sum(state_f)
                        if n_exc_f >= 1
                            for n in 1:N_oper_f
                                for k in 1:N_exp_term_f
                                    tmp = state_f[k + (n - 1) * N_exp_term_f]
                                    if tmp >= 1
                                        sum_ω += tmp * γ_list[n][k]
                                    end
                                end
                            end
                        end
                        row, col, val = pad_coo(- sum_ω * I_sup, N_he_tot, N_he_tot, idx + idx_f, idx + idx_f)
                        push!(localpart(L_row)[1], row...)
                        push!(localpart(L_col)[1], col...)
                        push!(localpart(L_val)[1], val...)
                    end
                    
                    # off-diagonal (boson)
                    state_neigh = copy(state_b)
                    for k in 1:N_exp_term_b
                        n_k = state_b[k]
                        if n_k >= 1
                            state_neigh[k] = n_k - 1
                            idx_neigh = ado2idx_b[state_neigh]
                            op = -1im * n_k * (c_list[k] * spreQ_b - conj(c_list[k]) * spostQ_b)
                            state_neigh[k] = n_k
                            for idx_f in 1:N_he_f
                                row, col, val = pad_coo(op, N_he_tot, N_he_tot, (idx + idx_f), (idx_neigh - 1) * N_he_f + idx_f)
                                push!(localpart(L_row)[1], row...)
                                push!(localpart(L_col)[1], col...)
                                push!(localpart(L_val)[1], val...)
                            end
                        end
                        if n_exc_b <= tier_b - 1
                            state_neigh[k] = n_k + 1
                            idx_neigh = ado2idx_b[state_neigh]
                            op = -1im * commQ_b
                            state_neigh[k] = n_k
                            for idx_f in 1:N_he_f
                                row, col, val = pad_coo(op, N_he_tot, N_he_tot, (idx + idx_f), (idx_neigh - 1) * N_he_f + idx_f)
                                push!(localpart(L_row)[1], row...)
                                push!(localpart(L_col)[1], col...)
                                push!(localpart(L_val)[1], val...)
                            end
                        end
                    end
                    if progressBar
                        put!(channel, true) # trigger a progress bar update
                    end
                    1 # Here, returning some number 1 and reducing it somehow (+) is necessary to make the distribution happen.
                end

                # fermion (n+1 & n-1 tier) superoperator
                @distributed (+) for idx_f in 1:N_he_f
                    state_f = idx2ado_f[idx_f]
                    n_exc_f = sum(state_f)
                    state_neigh = copy(state_f)
                    for n in 1:N_oper_f
                        for k in 1:(N_exp_term_f)
                            n_k = state_f[k + (n - 1) * N_exp_term_f]
                            if n_k >= 1
                                state_neigh[k + (n - 1) * N_exp_term_f] = n_k - 1
                                idx_neigh = ado2idx_f[state_neigh]
                                op = (-1) ^ eval(parity) * η_list[n][k] * spreQ_f[n] - (-1.0) ^ (n_exc_f - 1) * conj(η_list[(n % 2 == 0) ? (n-1) : (n+1)][k]) * spostQ_f[n]

                            elseif n_exc_f <= tier_f - 1
                                state_neigh[k + (n - 1) * N_exp_term_f] = n_k + 1
                                idx_neigh = ado2idx_f[state_neigh]
                                op = (-1) ^ eval(parity) * spreQd_f[n] + (-1.0) ^ (n_exc_f + 1) * spostQd_f[n]

                            else
                                continue
                            end
                            tmp_exc = sum(state_neigh[1:(k + (n - 1) * N_exp_term_f - 1)])
                            op *= -1im * (-1) ^ (tmp_exc)
                            state_neigh[k + (n - 1) * N_exp_term_f] = n_k

                            for idx_b in 1:N_he_b
                                idx = (idx_b - 1) * N_he_f
                                row, col, val = pad_coo(op, N_he_tot, N_he_tot, idx + idx_f, idx + idx_neigh)
                                push!(localpart(L_row)[1], row...)
                                push!(localpart(L_col)[1], col...)
                                push!(localpart(L_val)[1], val...)
                            end
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
        println("Constructing matrix...")
        L_he = sparse(vcat(L_row...), vcat(L_col...), vcat(L_val...), N_he_tot * sup_dim, N_he_tot * sup_dim)

        # add the liouville of system Hamiltonian term
        L_he += kron(sparse(I, N_he_tot, N_he_tot), -1im * (spre(Hsys) - spost(Hsys)))

        println("[DONE]")
        return new(L_he, tier_b, tier_f, Nsys, N_he_tot, N_he_b, N_he_f, sup_dim, parity, ado2idx_b_ordered, ado2idx_f_ordered)
    end
end
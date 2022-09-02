"""
# `M_Fermion <: AbstractHEOMMatrix`
Heom matrix for fermionic bath

## Fields
- `data::SparseMatrixCSC{ComplexF64, Int64}` : the sparse matrix
- `tier::Int`   : the tier (cutoff) for the bath
- `dim::Int`  : the dimension of system
- `N::Int`   : the number of total states
- `Nb::Int`   : the number of bosonic states (should be zero)
- `Nf::Int`   : the number of fermionic states
- `sup_dim::Int`: the dimension of system superoperator
- `parity::Symbol`: the parity of the density matrix
- `ado2idx::OrderedDict{Vector{Int}, Int}`: the ADO-to-index dictionary

## Constructor
`M_Fermion(Hsys, tier, bath, parity; [progressBar])`

- `Hsys::AbstractMatrix` : The system Hamiltonian
- `tier::Int` : the tier (cutoff) for the bath
- `bath::FermionicBath` : an object for the fermionic bath correlation
- `parity::Symbol` : The parity symbol of the density matrix (either `:odd` or `:even`). Defaults to `:even`.
- `progressBar::Bool` : Display progress bar during the process or not. Defaults to `true`.
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
    const ado2idx::OrderedDict{Vector{Int}, Int}
    
    function M_Fermion(        
            Hsys::AbstractMatrix,
            tier::Int,
            bath::FermionicBath,
            parity::Symbol=:even;
            progressBar::Bool=true
        )

        if (parity != :even) || (parity != :odd)
            error("The parity symbol of density matrix should be either \":odd\" or \":even\".")
        end

        Nsys,   = size(Hsys)
        sup_dim = Nsys ^ 2
        I_sup   = sparse(I, sup_dim, sup_dim)

        η_list = bath.η_list
        γ_list = bath.γ_list
        Coup_Ops   = bath.coupOP
        N_oper     = bath.N_oper
        N_exp_term = bath.N_term
            
        dims    = [2 for i in 1:(N_exp_term * N_oper)]
    
        spreQ   = spre.(Coup_Ops)
        spostQ  = spost.(Coup_Ops)
        spreQd  = spre.(adjoint.(Coup_Ops))
        spostQd = spost.(adjoint.(Coup_Ops))

        # get ADOs dictionary
        N_he, ado2idx_ordered, idx2ado = ADOs_dictionary(dims, tier)
        ado2idx = Dict(ado2idx_ordered)

        # start to construct the matrix
        L_row = distribute([Int[] for _ in procs()])
        L_col = distribute([Int[] for _ in procs()])
        L_val = distribute([ComplexF64[] for _ in procs()])
        channel = RemoteChannel(() -> Channel{Bool}(), 1) # for updating the progress bar

        println("Start constructing hierarchy matrix...(using $(nprocs()) processors)")
        if progressBar
            prog = Progress(N_he; desc="Processing: ", PROGBAR_OPTIONS...)
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
                @distributed (+) for idx in 1:N_he
                    state = idx2ado[idx]
                    n_exc = sum(state)
                    sum_ω = 0.0
                    if n_exc >= 1
                        for n in 1:N_oper
                            for k in 1:(N_exp_term)
                                tmp = state[k + (n - 1) * N_exp_term]
                                if tmp >= 1
                                    sum_ω += tmp * γ_list[n][k]
                                end
                            end
                        end
                    end
                    row, col, val = pad_coo(- sum_ω * I_sup, N_he, N_he, idx, idx)
                    push!(localpart(L_row)[1], row...)
                    push!(localpart(L_col)[1], col...)
                    push!(localpart(L_val)[1], val...)

                    state_neigh = copy(state)
                    for n in 1:N_oper
                        for k in 1:(N_exp_term)
                            n_k = state[k + (n - 1) * N_exp_term]
                            if n_k >= 1
                                state_neigh[k + (n - 1) * N_exp_term] = n_k - 1
                                idx_neigh = ado2idx[state_neigh]
                                op = (-1) ^ eval(parity) * η_list[n][k] * spreQ[n] - (-1.0) ^ (n_exc - 1) * conj(η_list[(n % 2 == 0) ? (n-1) : (n+1)][k]) * spostQ[n]

                            elseif n_exc <= tier - 1
                                state_neigh[k + (n - 1) * N_exp_term] = n_k + 1
                                idx_neigh = ado2idx[state_neigh]
                                op = (-1) ^ eval(parity) * spreQd[n] + (-1.0) ^ (n_exc + 1) * spostQd[n]

                            else
                                continue
                            end

                            tmp_exc = sum(state_neigh[1:(k + (n - 1) * N_exp_term - 1)])
                            row, col, val = pad_coo(-1im * (-1) ^ (tmp_exc) * op, N_he, N_he, idx, idx_neigh)
                            push!(localpart(L_row)[1], row...)
                            push!(localpart(L_col)[1], col...)
                            push!(localpart(L_val)[1], val...)

                            state_neigh[k + (n - 1) * N_exp_term] = n_k
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
        L_he = sparse(vcat(L_row...), vcat(L_col...), vcat(L_val...), N_he * sup_dim, N_he * sup_dim)

        # add the liouville of system Hamiltonian term
        L_he += kron(sparse(I, N_he_tot, N_he_tot), -1im * (spre(Hsys) - spost(Hsys)))
        
        println("[DONE]")
        return new(L_he, tier, Nsys, N_he, 0, N_he, sup_dim, parity, ado2idx_ordered)
    end
end
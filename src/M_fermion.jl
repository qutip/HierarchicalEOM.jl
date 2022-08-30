"""
# `M_Fermion <: AbstractHEOMMatrix`
Heom matrix for fermionic bath

## Fields
- `data::SparseMatrixCSC{ComplexF64, Int64}` : the sparse matrix
- `tier::Int`   : the tier (cutoff) for the bath
- `N_sys::Int`  : the dimension of system
- `N_he::Int`   : the number of states
- `sup_dim::Int`: the dimension of system superoperator
- `ado2idx::OrderedDict{Vector{Int}, Int}`: the ADO-to-index dictionary

## Constructor
`M_Fermion(Hsys, tier, η_list, γ_list, Coup_Ops; [spectral, progressBar])`

- `Hsys::AbstractMatrix` : The system Hamiltonian
- `tier::Int` : the tier (cutoff) for the bath
- `η_list::Vector{Vector{Tv<:Number}}` : the coefficient ``\\eta_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `γ_list::Vector{Vector{Ti<:Number}}` : the coefficient ``\\gamma_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `Coup_Ops::Vector` : Operator list describing the coupling between system and bath.
- `spectral::Bool` : Decide whether to calculate spectral density or not. Defaults to `false`.
- `progressBar::Bool` : Display progress bar during the process or not. Defaults to `true`.
"""
mutable struct M_Fermion <: AbstractHEOMMatrix
    data::SparseMatrixCSC{ComplexF64, Int64}
    tier::Int
    N_sys::Int
    N_he::Int
    sup_dim::Int
    ado2idx::OrderedDict{Vector{Int}, Int}
    
    function M_Fermion(        
            Hsys::AbstractMatrix,
            tier::Int,
            η_list::Vector{Vector{Tv}},
            γ_list::Vector{Vector{Ti}},
            Coup_Ops::Vector;
            spectral::Bool=false,
            progressBar::Bool=true
        ) where {Ti,Tv <: Number}

        # check if the length of η_list, γ_list, and Coup_Ops are valid
        N_bath  = length(Coup_Ops)
        N_exp_term = length(η_list[1])
        if (N_bath != length(η_list)) || (N_bath != length(γ_list))
            error("The length of \'η_list\', \'γ_list\', and \'Coup_Ops\' should all be the same.")
        end
        for i in 1:N_bath
            if N_exp_term != length(η_list[i]) 
                error("The length of each vector in \'η_list\' are wrong.")
            end
            if N_exp_term != length(γ_list[i]) 
                error("The length of each vector in \'γ_list\' are wrong.")
            end
        end
            
        Nsys,   = size(Hsys)
        dims    = [2 for i in 1:(N_exp_term * N_bath)]
        sup_dim = Nsys ^ 2
        I_sup   = sparse(I, sup_dim, sup_dim)
    
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
                        for n in 1:N_bath
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
                    for n in 1:N_bath
                        for k in 1:(N_exp_term)
                            n_k = state[k + (n - 1) * N_exp_term]
                            if n_k >= 1
                                state_neigh[k + (n - 1) * N_exp_term] = n_k - 1
                                idx_neigh = ado2idx[state_neigh]
                                op = (-1) ^ spectral * η_list[n][k] * spreQ[n] - (-1.0) ^ (n_exc - 1) * conj(η_list[(n % 2 == 0) ? (n-1) : (n+1)][k]) * spostQ[n]

                            elseif n_exc <= tier - 1
                                state_neigh[k + (n - 1) * N_exp_term] = n_k + 1
                                idx_neigh = ado2idx[state_neigh]
                                op = (-1) ^ (spectral) * spreQd[n] + (-1.0) ^ (n_exc + 1) * spostQd[n]

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

        # add the free Hamiltonian evolution term
        L_he += kron(sparse(I, N_he_tot, N_he_tot), -1im * (spre(Hsys) - spost(Hsys)))
        
        println("[DONE]")
        return new(L_he, tier, Nsys, N_he, sup_dim, ado2idx_ordered)
    end
end
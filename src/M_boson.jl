"""
# `M_boson <: AbstractHEOMMatrix`
Heom matrix for bosonic bath

## Fields
- `data::SparseMatrixCSC{ComplexF64, Int64}` : the sparse matrix
- `tier::Int`   : the tier (cutoff) for the bath
- `N_sys::Int`  : the dimension of system
- `N_he::Int`   : the number of states
- `sup_dim::Int`: the dimension of system superoperator
- `ADOs::OrderedDict{Vector{Int}, Int}`: the ADOs dictionary

## Constructor
`M_boson(Hsys, tier, η_list, γ_list, Coup_Op; [Jump_Ops, liouville])`

- `Hsys::AbstractMatrix` : The system Hamiltonian
- `tier::Int` : the tier (cutoff) for the bath
- `η_list::Vector{Tv<:Number}` : the coefficient ``\\eta_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `γ_list::Vector{Ti<:Number}` : the coefficient ``\\gamma_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `Coup_Op::AbstractMatrix` : Operator describing the coupling between system and bath.
- `Jump_Ops::Vector` : The collapse (jump) operators to add when calculating liouvillian in lindblad term (only if `liouville=true`). Defaults to empty vector `[]`.
- `liouville::Bool` : Add liouvillian to the matrix or not. Defaults to `true`.
- `progressBar::Bool` : Display progress bar during the process or not. Defaults to `true`.
"""
mutable struct M_boson <: AbstractHEOMMatrix
    data::SparseMatrixCSC{ComplexF64, Int64}
    tier::Int
    N_sys::Int
    N_he::Int
    sup_dim::Int
    ADOs::OrderedDict{Vector{Int}, Int}
    
    function M_boson(        
            Hsys::AbstractMatrix,
            tier::Int,
            η_list::Vector{Tv},
            γ_list::Vector{Ti},
            Coup_Op::AbstractMatrix;
            Jump_Ops::Vector=[],   # only when liouville is set to true
            liouville::Bool=true,
            progressBar::Bool=true
        ) where {Ti,Tv <: Number}

        # check if the length of η_list and γ_list are valid
        N_exp_term = length(η_list)
        if N_exp_term != length(γ_list)
            error("The length of \'η_list\' and \'γ_list\' should be the same.")
        end
    
        Nsys,   = size(Hsys)
        dims    = [(tier + 1) for i in 1:N_exp_term]
        sup_dim = Nsys ^ 2
        I_sup   = sparse(I, sup_dim, sup_dim)

        spreQ  = spre(Coup_Op)
        spostQ = spost(Coup_Op)
        commQ  = spreQ - spostQ

        # get ADOs dictionary
        N_he, he2idx_ordered, idx2he = ADOs_dictionary(dims, tier)
        he2idx = Dict(he2idx_ordered)

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
                    state = idx2he[idx]
                    n_exc = sum(state)
                    sum_ω = 0.0
                    if n_exc >= 1
                        for k in 1:N_exp_term
                            if state[k] > 0
                                sum_ω += state[k] * γ_list[k]
                            end
                        end
                    end
                    row, col, val = pad_coo(- sum_ω * I_sup, N_he, N_he, idx, idx)
                    push!(localpart(L_row)[1], row...)
                    push!(localpart(L_col)[1], col...)
                    push!(localpart(L_val)[1], val...)

                    state_neigh = copy(state)
                    for k in 1:N_exp_term
                        n_k = state[k]
                        if n_k >= 1
                            state_neigh[k] = n_k - 1
                            idx_neigh = he2idx[state_neigh]
                            op = -1im * n_k * (η_list[k] * spreQ - conj(η_list[k]) * spostQ)
                            row, col, val = pad_coo(op, N_he, N_he, idx, idx_neigh)
                            push!(localpart(L_row)[1], row...)
                            push!(localpart(L_col)[1], col...)
                            push!(localpart(L_val)[1], val...)
                            
                            state_neigh[k] = n_k
                        end
                        if n_exc <= tier - 1
                            state_neigh[k] = n_k + 1
                            idx_neigh = he2idx[state_neigh]
                            op = -1im * commQ
                            row, col, val = pad_coo(op, N_he, N_he, idx, idx_neigh)
                            push!(localpart(L_row)[1], row...)
                            push!(localpart(L_col)[1], col...)
                            push!(localpart(L_val)[1], val...)
                            
                            state_neigh[k] = n_k
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

        if liouville
            println("Adding liouvillian...")
            L_he += kron(sparse(I, N_he, N_he), liouvillian(Hsys, Jump_Ops, progressBar))
        end
        
        println("[DONE]")
        return new(L_he, tier, Nsys, N_he, sup_dim, he2idx_ordered)
    end
end
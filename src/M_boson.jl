"""
# `M_Boson <: AbstractHEOMMatrix`
Heom matrix for bosonic bath

## Fields
- `data::SparseMatrixCSC{ComplexF64, Int64}` : the sparse matrix
- `tier::Int` : the tier (cutoff) for the bath
- `dim::Int`  : the dimension of system
- `N::Int`  : the number of total states
- `Nb::Int` : the number of bosonic states
- `Nf::Int` : the number of fermionic states (should be zero)
- `sup_dim::Int` : the dimension of system superoperator
- `parity::Symbol` : the parity of the density matrix (restrict to `:none` for boson)
- `ado2idx::OrderedDict{Vector{Int}, Int}` : the ADO-to-index dictionary

## Constructor
`M_Boson(Hsys, tier, bath; [progressBar])`

- `Hsys::AbstractMatrix` : The system Hamiltonian
- `tier::Int` : the tier (cutoff) for the bath
- `bath::BosonicBath` : an object for the bosonic bath correlation
- `progressBar::Bool` : Display progress bar during the process or not. Defaults to `true`.
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
    const ado2idx::OrderedDict{Vector{Int}, Int}
    
    function M_Boson(        
            Hsys::AbstractMatrix,
            tier::Int,
            bath::BosonicBath;
            progressBar::Bool=true
        )

        c_list = bath.η_list
        ν_list = bath.γ_list
        Coup_Op = bath.coupOP
        N_exp_term = bath.N_term
    
        Nsys,   = size(Hsys)
        dims    = [(tier + 1) for i in 1:N_exp_term]
        sup_dim = Nsys ^ 2
        I_sup   = sparse(I, sup_dim, sup_dim)

        spreQ  = spre(Coup_Op)
        spostQ = spost(Coup_Op)
        commQ  = spreQ - spostQ

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
                        for k in 1:N_exp_term
                            if state[k] > 0
                                sum_ω += state[k] * ν_list[k]
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
                            idx_neigh = ado2idx[state_neigh]
                            op = -1im * n_k * (c_list[k] * spreQ - conj(c_list[k]) * spostQ)
                            row, col, val = pad_coo(op, N_he, N_he, idx, idx_neigh)
                            push!(localpart(L_row)[1], row...)
                            push!(localpart(L_col)[1], col...)
                            push!(localpart(L_val)[1], val...)
                            
                            state_neigh[k] = n_k
                        end
                        if n_exc <= tier - 1
                            state_neigh[k] = n_k + 1
                            idx_neigh = ado2idx[state_neigh]
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

        # add the liouville of system Hamiltonian term
        L_he += kron(sparse(I, N_he, N_he), -1im * (spre(Hsys) - spost(Hsys)))
        
        println("[DONE]")
        return new(L_he, tier, Nsys, N_he, N_he, 0, sup_dim, :none, ado2idx_ordered)
    end
end
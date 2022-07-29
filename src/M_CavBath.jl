"""
# `M_CavBath <: AbstractHEOMMatrix`
Heom matrix with setting the single mode cavity as bosonic bath

## Fields
- `data::SparseMatrixCSC{ComplexF64, Int64}` : the sparse matrix
- `tier_b::Int`   : the tier (cutoff) for bosonic bath
- `tier_f::Int`   : the tier (cutoff) for fermionic bath
- `N_sys::Int`  : the dimension of system
- `N_he::Int`   : the number of total states
- `N_he_b::Int`   : the number of bosonic states
- `N_he_f::Int`   : the number of fermionic states
- `sup_dim::Int`: the dimension of system superoperator
- `ados_b::OrderedDict{Vector{Int}, Int}`: the bosonic ados dictionary
- `ados_f::OrderedDict{Vector{Int}, Int}`: the fermionic ados dictionary

## Constructor
`M_CavBath(Hsys, tier_b, tier_f, c_list, ν_list, η_list, γ_list, Coup_Op_b, Coup_Op_f; [Jump_Ops, spectral, liouville])`

- `Hsys::Union{AbstractMatrix, AbstractOperator}` : The system Hamiltonian
- `tier_b::Int` : the tier (cutoff) for the bosonic bath
- `tier_f::Int` : the tier (cutoff) for the fermionic bath
- `c_list::Vector{Ti<:Number}` : the coefficient ``c_i`` in bosonic bath correlation functions (``\\sum_i c_i e^{-\\nu_i t}``).
- `ν_list::Vector{Tj<:Number}` : the coefficient ``\\nu_i`` in bosonic bath correlation functions (``\\sum_i c_i e^{-\\nu_i t}``)
- `η_list::Vector{Vector{Tk<:Number}}` : the coefficient ``\\eta_i`` in fermionic bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `γ_list::Vector{Vector{Tl<:Number}}` : the coefficient ``\\gamma_i`` in fermionic bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
- `Coup_Op_b::Vector` : Operator list describing the coupling between system and bosonic bath.
- `Coup_Op_f::Vector` : Operator list describing the coupling between system and fermionic bath.
- `Jump_Ops::Vector` : The collapse (jump) operators to add when calculating liouvillian in lindblad term (only if `liouville=true`). Defaults to empty vector `[]`.
- `spectral::Bool` : Decide whether to calculate spectral density or not. Defaults to `false`.
- `liouville::Bool` : Add liouvillian to the matrix or not. Defaults to `true`.
"""
mutable struct M_CavBath <: AbstractHEOMMatrix
    data::SparseMatrixCSC{ComplexF64, Int64}
    tier_b::Int
    tier_f::Int
    N_sys::Int
    N_he::Int
    N_he_b::Int
    N_he_f::Int
    sup_dim::Int
    ados_b::OrderedDict{Vector{Int}, Int}
    ados_f::OrderedDict{Vector{Int}, Int}
    
    function M_CavBath(        
            Hsys::Union{AbstractMatrix, AbstractOperator},
            tier_b::Int,
            tier_f::Int,
            c_list::Vector{Ti},
            ν_list::Vector{Tj},
            η_list::Vector{Vector{Tk}},
            γ_list::Vector{Vector{Tl}},
            Coup_Op_b::Vector,
            Coup_Op_f::Vector;
            Jump_Ops::Vector=[],  # only when liouville is set to true
            spectral::Bool=false,
            liouville::Bool=true
        ) where {Ti,Tj,Tk,Tl <: Number}

        # check if the length of c_list, ν_list, and Coup_Ops are valid
        N_exp_term_b = length(Coup_Op_b)
        if (N_exp_term_b != length(c_list)) || (N_exp_term_b != length(ν_list))
            error("The length of \'Coup_Op_b\', \'c_list\', and \'ν_list\' should be the same.")
        end

        # check if the length of η_list, γ_list, and Coup_Ops are valid
        N_bath_f  = length(Coup_Op_f)
        N_exp_term_f = length(η_list[1])
        if (N_bath_f != length(η_list)) || (N_bath_f != length(γ_list))
            error("The length of \'η_list\', \'γ_list\', and \'f_Coup_Op\' should all be the same.")
        end
        for i in 1:N_bath_f
            if N_exp_term_f != length(η_list[i]) 
                error("The length of each vector in \'η_list\' are wrong.")
            end
            if N_exp_term_f != length(γ_list[i]) 
                error("The length of each vector in \'γ_list\' are wrong.")
            end
        end
            
        Nsys,   = size(Hsys)
        sup_dim = Nsys ^ 2
        I_sup   = sparse(I, sup_dim, sup_dim)

        dims_b    = [(tier_b + 1) for i in 1:N_exp_term_b]
        dims_f    = [2 for i in 1:(N_exp_term_f * N_bath_f)]
    
        spreQ_b   = spre.(Coup_Op_b)
        spostQ_b  = spost.(Coup_Op_b)
        commQ_b   = spreQ_b .- spostQ_b
        spreQ_f   = spre.(Coup_Op_f)
        spostQ_f  = spost.(Coup_Op_f)
        spreQd_f  = spre.(dagger.(Coup_Op_f))
        spostQd_f = spost.(dagger.(Coup_Op_f))

        # get Ados dictionary
        N_he_b, he2idx_b_ordered, idx2he_b = Ados_dictionary(dims_b, tier_b)
        N_he_f, he2idx_f_ordered, idx2he_f = Ados_dictionary(dims_f, tier_f)
        N_he_tot = N_he_b * N_he_f
        he2idx_b = Dict(he2idx_b_ordered)
        he2idx_f = Dict(he2idx_f_ordered)

        # start to construct the matrix L_he
        println("Start constructing process...")
        L_he = spzeros(ComplexF64, N_he_tot * sup_dim, N_he_tot * sup_dim)

        for idx_b in 1:N_he_b
            # diagonal (boson)
            sum_ω   = 0.0
            state_b = idx2he_b[idx_b]
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
                state_f = idx2he_f[idx_f]
                n_exc_f = sum(state_f)
                if n_exc_f >= 1
                    for n in 1:N_bath_f
                        for k in 1:N_exp_term_f
                            tmp = state_f[k + (n - 1) * N_exp_term_f]
                            if tmp >= 1
                                sum_ω += tmp * γ_list[n][k]
                            end
                        end
                    end
                end
                L_he += pad_csc(- sum_ω * I_sup, N_he_tot, N_he_tot, idx + idx_f, idx + idx_f)
            end
            
            # off-diagonal (boson)
            state_neigh = copy(state_b)
            for k in 1:N_exp_term_b
                n_k = state_b[k]
                if n_k >= 1
                    state_neigh[k] = n_k - 1
                    idx_neigh = he2idx_b[state_neigh]
                    op = c_list[k] * spreQ_b[k] - c_list[(k % 2 == 0) ? (k-1) : (k+1)] * spostQ_b[k]
                    state_neigh[k] = n_k
                    for idx_f in 1:N_he_f
                        L_he += pad_csc(op, N_he_tot, N_he_tot, (idx + idx_f), (idx_neigh - 1) * N_he_f + idx_f)
                    end
                end
                if n_exc_b <= tier_b - 1
                    state_neigh[k] = n_k + 1
                    idx_neigh = he2idx_b[state_neigh]
                    op = -1im * commQ_b[k]
                    state_neigh[k] = n_k
                    for idx_f in 1:N_he_f
                        L_he += pad_csc(op, N_he_tot, N_he_tot, (idx + idx_f), (idx_neigh - 1) * N_he_f + idx_f)
                    end
                end
            end
        end

        # fermion (n+1 & n-1 tier) superoperator
        for idx_f in 1:N_he_f
            state_f = idx2he_f[idx_f]
            n_exc_f = sum(state_f)
            state_neigh = copy(state_f)
            for n in 1:N_bath_f
                for k in 1:(N_exp_term_f)
                    n_k = state_f[k + (n - 1) * N_exp_term_f]
                    if n_k >= 1
                        state_neigh[k + (n - 1) * N_exp_term_f] = n_k - 1
                        idx_neigh = he2idx_f[state_neigh]
                        op = (-1) ^ spectral * η_list[n][k] * spreQ_f[n] - (-1.0) ^ (n_exc_f - 1) * conj(η_list[(n % 2 == 0) ? (n-1) : (n+1)][k]) * spostQ_f[n]

                    elseif n_exc_f <= tier_f - 1
                        state_neigh[k + (n - 1) * N_exp_term_f] = n_k + 1
                        idx_neigh = he2idx_f[state_neigh]
                        op = (-1) ^ (spectral) * spreQd_f[n] + (-1.0) ^ (n_exc_f + 1) * spostQd_f[n]

                    else
                        continue
                    end
                    tmp_exc = sum(state_neigh[1:(k + (n - 1) * N_exp_term_f - 1)])
                    op *= -1im * (-1) ^ (tmp_exc)
                    state_neigh[k + (n - 1) * N_exp_term_f] = n_k

                    for idx_b in 1:N_he_b
                        idx = (idx_b - 1) * N_he_f
                        L_he += pad_csc(op, N_he_tot, N_he_tot, idx + idx_f, idx + idx_neigh)
                    end
                end
            end
        end

        if liouville
            print("Construct Liouvillian...")
            L_he += kron(sparse(I, N_he_tot, N_he_tot), liouvillian(Hsys, Jump_Ops))
        end

        println("[DONE]")
        return new(L_he, tier_b, tier_f, Nsys, N_he_tot, N_he_b, N_he_f, sup_dim, he2idx_b_ordered, he2idx_f_ordered)
    end
end
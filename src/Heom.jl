module Heom
    import Base: size
    import LinearAlgebra: eigvals, I
    import DifferentialEquations: ODEProblem, init, Tsit5, step!
    import SparseArrays: sparse, spzeros, SparseMatrixCSC, SparseVector, AbstractSparseMatrix
    import QuantumOptics: AbstractOperator, dagger
    import ProgressMeter: Progress, next!

    export AbstractHEOMMatrix, M_fermion, M_boson, evolution, pade_NmN, Correlation, spre, spost, liouvillian

    PROGBAR_OPTIONS = Dict(:barlen=>20, :color=>:green, :showspeed=>true)

    fb(m) = 2 * m + 1

    spre(q::AbstractMatrix)        = kron(Matrix(I, size(q)[1], size(q)[1]), q)
    spre(q::AbstractOperator)      = sparse(kron(sparse(I, size(q)[1], size(q)[1]), q.data))
    spre(q::AbstractSparseMatrix)  = sparse(kron(sparse(I, size(q)[1], size(q)[1]), q))
    spost(q::AbstractMatrix)       = kron(transpose(q), Matrix(I, size(q)[1], size(q)[1]))
    spost(q::AbstractOperator)     = sparse(kron(transpose(q.data), sparse(I, size(q)[1], size(q)[1])))
    spost(q::AbstractSparseMatrix) = sparse(kron(transpose(q), sparse(I, size(q)[1], size(q)[1])))

    function δ(j, k)
        if j == k
            return 1
        else
            return 0
        end
    end

    function i_eigval(r, p, s)
        pf = []
        A_matrix = zeros(r, r)
        for j in 1:r
            for k in 1:r
                A_matrix[j, k] = (δ(j, k+1) + δ(j, k-1)) / √((fb(j-1) + p) * (fb(k-1) + p))
            end
        end

        for val in eigvals(A_matrix)[1:s]
            push!(pf, -2 / val)
        end

        return pf
    end

    #thoss's spectral (N-1/N) pade
    function pade_NmN(lmax,fermion::Bool=true)

        if fermion
            local ϵ = i_eigval(2 * lmax    , 0, lmax)
            local χ = i_eigval(2 * lmax - 1, 2, lmax - 1)
            prefactor = 0.5 * lmax * fb(lmax)
        else
            local ϵ = i_eigval(2 * lmax    , 2, lmax)
            local χ = i_eigval(2 * lmax - 1, 4, lmax - 1)
            prefactor = 0.5 * lmax * (fb(lmax) + 2)
        end

        η_list = []
        for j in 1:lmax
            term = prefactor
            for k1 in 1:(lmax - 1)
                term *= (χ[k1] ^ 2 - ϵ[j] ^ 2) / (ϵ[k1] ^ 2 - ϵ[j] ^ 2 + δ(j, k1))
            end

            term /= (ϵ[lmax] ^ 2 - ϵ[j] ^ 2 + δ(j, lmax))

            push!(η_list, term)
        end

        return append!([0.0], η_list), append!([0.0], ϵ)
    end

    function f_approx(x, η, ϵ, lmax)
        f = 0.5
        for l in 2:(lmax + 1)
            f -= 2 * η[l] * x / (x ^ 2 + ϵ[l] ^ 2)
        end
        return f
    end 

    # Correlation function
    function Correlation(σ, μ, η, ϵ, β, W, Γ, lmax)
        local η_list = [0.5 * Γ * W * f_approx(1.0im * β * W, η, ϵ, lmax)]
        local γ_list = [W - σ * 1.0im * μ]

        if lmax > 0
            for l in 2:(lmax + 1)
                append!(η_list, -1.0im * (η[l] / β) * Γ * W ^ 2 / (-(ϵ[l] / β) ^ 2 + W ^ 2))
                append!(γ_list, ϵ[l] / β - σ * 1.0im * μ)
            end
        end
        return η_list, γ_list
    end

    # generate index to state vector
    function state_number(dims::Vector{Int}, N_exc::Int)
        len = length(dims)
        state = zeros(Int, len)
        result = [copy(state)]
        nexc = 0

        while true
            idx = len
            state[end] += 1
            nexc += 1
            if state[idx] < dims[idx]
                push!(result, copy(state))
            end
            while (nexc == N_exc) || (state[idx] == dims[idx])
                #state[idx] = 0
                idx -= 1
                if idx < 1
                    return result
                end

                nexc -= state[idx + 1] - 1
                state[idx + 1] = 0
                state[idx] += 1
                if state[idx] < dims[idx]
                    push!(result, copy(state))
                end
            end
        end
    end

    function Ados_dictionary(dims::Vector{Int}, N_exc::Int)
        state2idx = Dict{Vector{Int}, Int}()
        idx2state = state_number(dims, N_exc)
        for (idx, state) in enumerate(idx2state)
            state2idx[state] = idx
        end

        return length(idx2state), state2idx, idx2state
    end

    function csr2csc(A::AbstractSparseMatrix)
        n = A.n
        ptr = A.colptr

        # return row_idx, col_idx, value
        return A.rowval, [i for i in 1:n for j in ptr[i]:(ptr[i+1]-1)], A.nzval
    end

    function pad_csc(A::SparseMatrixCSC{T, Int64}, row_scale::Int, col_scale::Int, row_idx=1::Int, col_idx=1::Int) where {T<:Number}
        (M, N) = size(A)

        # deal with values
        values = A.nzval
        if length(values) == 0
            return sparse([M * row_scale], [N * col_scale], [0.0im])
        else
            if T != ComplexF64
                values = convert.(ComplexF64, values)
            end

            # deal with colptr
            local ptrLen::Int         = N * col_scale + 1
            local ptrIn::Vector{Int}  = A.colptr
            local ptrOut::Vector{Int} = fill(1, ptrLen)
            if col_idx == 1
                ptrOut[1:(N+1)]   .= ptrIn            
                ptrOut[(N+2):end] .= ptrIn[end]

            elseif col_idx == col_scale         
                ptrOut[(ptrLen-N):end] .= ptrIn

            elseif (col_idx < col_scale) && (col_idx > 1)
                tmp1 = (col_idx - 1) * N + 1
                tmp2 = tmp1 + N
                ptrOut[tmp1:tmp2] .= ptrIn
                ptrOut[(tmp2+1):end] .= ptrIn[end]

            else
                error("col_idx must be \'>= 1\' and \'<= col_scale\'")
            end

            # deal with rowval
            if (row_idx > row_scale) || (row_idx < 1)
                error("row_idx must be \'>= 1\' and \'<= row_scale\'")
            end
            tmp1 = (row_idx - 1) * N

            return SparseMatrixCSC(
                M * row_scale,
                N * col_scale,
                ptrOut,
                A.rowval .+ tmp1, 
                values,
            )
        end
    end

    function liouvillian(Hsys, Jump_Ops::Vector=[], progressBar::Bool=true)
        
        N, = size(Hsys)

        L = -1im * (spre(Hsys) - spost(Hsys))
        if progressBar
            prog = Progress(length(Jump_Ops) + 1, start=1; desc="Construct Liouvillian     : ", PROGBAR_OPTIONS...)
        end
        for J in Jump_Ops
            L += spre(J) * spost(J') - 0.5 * (spre(J' * J) + spost(J' * J))
            if progressBar
                next!(prog)
            end
        end
        return L
    end

    abstract type AbstractHEOMMatrix end
    size(A::AbstractHEOMMatrix) = size(A.data)

    """
    # `M_fermion <: AbstractHEOMMatrix`
    Heom matrix for fermionic bath

    ## Fields
    - `data::SparseMatrixCSC{ComplexF64, Int64}` : the sparse matrix
    - `tier::Int`   : the tier (cutoff) for the bath
    - `N_sys::Int`  : the dimension of system
    - `N_he::Int`   : the number of states
    - `sup_dim::Int`: the dimension of system superoperator

    ## Constructor
    `M_fermion(Hsys, tier, η_list, γ_list, Coup_Ops; [Jump_Ops, spectral, liouville, progressBar])`

    - `Hsys::Union{AbstractMatrix, AbstractOperator}` : The system Hamiltonian
    - `tier::Int` : the tier (cutoff) for the bath
    - `η_list::Vector{Vector{Tv<:Number}}` : the coefficient ``\\eta_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
    - `γ_list::Vector{Vector{Ti<:Number}}` : the coefficient ``\\gamma_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
    - `Coup_Ops::Vector` : Operator list describing the coupling between system and bath.
    - `Jump_Ops::Vector` : The collapse (jump) operators to add when calculating liouvillian in lindblad term (only if `liouville=true`). Defaults to empty vector `[]`.
    - `spectral::Bool` : Decide whether to calculate spectral density or not. Defaults to `false`.
    - `liouville::Bool` : Add liouvillian to the matrix or not. Defaults to `true`.
    - `progressBar::Bool` : Display progress bar during the process or not. Defaults to `true`.
    """
    mutable struct M_fermion <: AbstractHEOMMatrix
        data::SparseMatrixCSC{ComplexF64, Int64}
        tier::Int
        N_sys::Int
        N_he::Int
        sup_dim::Int
        
        function M_fermion(        
                Hsys::Union{AbstractMatrix, AbstractOperator},
                tier::Int,
                η_list::Vector{Vector{Tv}},
                γ_list::Vector{Vector{Ti}},
                Coup_Ops::Vector;
                Jump_Ops::Vector=[],  # only when liouville is set to true
                spectral::Bool=false,
                liouville::Bool=true,
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
            spreQd  = spre.(dagger.(Coup_Ops))
            spostQd = spost.(dagger.(Coup_Ops))

            # get Ados dictionary
            N_he, he2idx, idx2he = Ados_dictionary(dims, tier)

            # start to construct the matrix L_he
            print("Start constructing process...")
            L_he = spzeros(ComplexF64, N_he * sup_dim, N_he * sup_dim)

            if liouville || progressBar
                print("\n")
            end
            flush(stdout)
            if progressBar
                prog = Progress(N_he; desc="Construct hierarchy matrix: ", PROGBAR_OPTIONS...)
            end
            for idx in 1:N_he
                state = idx2he[idx]
                n_exc = sum(state)
                sum_ω = 0.0
                for n in 1:N_bath
                    for k in 1:(N_exp_term)
                        if n_exc >= 1
                            sum_ω += state[k + (n - 1) * N_exp_term] * γ_list[n][k]
                        end
                    end
                end
                L_he += pad_csc(- sum_ω * I_sup, N_he, N_he, idx, idx)

                state_neigh = copy(state)
                for n in 1:N_bath
                    for k in 1:(N_exp_term)
                        n_k = state[k + (n - 1) * N_exp_term]
                        if n_k >= 1
                            state_neigh[k + (n - 1) * N_exp_term] = n_k - 1
                            idx_neigh = he2idx[state_neigh]
                            op = (-1) ^ spectral * η_list[n][k] * spreQ[n] - (-1.0) ^ (n_exc - 1) * conj(η_list[(n % 2 == 0) ? (n-1) : (n+1)][k]) * spostQ[n]

                        elseif n_exc <= tier - 1
                            state_neigh[k + (n - 1) * N_exp_term] = n_k + 1
                            idx_neigh = he2idx[state_neigh]
                            op = (-1) ^ (spectral) * spreQd[n] + (-1.0) ^ (n_exc + 1) * spostQd[n]

                        else
                            continue
                        end

                        tmp_exc = sum(state_neigh[1:(k + (n - 1) * N_exp_term - 1)])
                        L_he += pad_csc(-1im * (-1) ^ (tmp_exc) * op, N_he, N_he, idx, idx_neigh)
                        state_neigh[k + (n - 1) * N_exp_term] = n_k
                    end
                end
                if progressBar
                    next!(prog)
                end
            end

            if liouville
                printstyled("Construct Liouvillian...", color=:green)
                flush(stdout)
                L_he += kron(sparse(I, N_he, N_he), liouvillian(Hsys, Jump_Ops, progressBar))
            end
            
            println("[DONE]\n")
            flush(stdout)
            return new(L_he, tier, Nsys, N_he, sup_dim)
        end
    end

    """
    # `M_boson <: AbstractHEOMMatrix`
    Heom matrix for bosonic bath

    ## Fields
    - `data::SparseMatrixCSC{ComplexF64, Int64}` : the sparse matrix
    - `tier::Int`   : the tier (cutoff) for the bath
    - `N_sys::Int`  : the dimension of system
    - `N_he::Int`   : the number of states
    - `sup_dim::Int`: the dimension of system superoperator

    ## Constructor
    `M_boson(Hsys, tier, η_list, γ_list, Coup_Op; [Jump_Ops, liouville, progressBar])`

    - `Hsys::Union{AbstractMatrix, AbstractOperator}` : The system Hamiltonian
    - `tier::Int` : the tier (cutoff) for the bath
    - `η_list::Vector{Vector{Tv<:Number}}` : the coefficient ``\\eta_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
    - `γ_list::Vector{Vector{Ti<:Number}}` : the coefficient ``\\gamma_i`` in bath correlation functions (``\\sum_i \\eta_i e^{-\\gamma_i t}``).
    - `Coup_Op::Vector` : Operator describing the coupling between system and bath.
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
        
        function M_boson(        
                Hsys::Union{AbstractMatrix, AbstractOperator},
                tier::Int,
                η_list::Vector{Tv},
                γ_list::Vector{Ti},
                Coup_Op;
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

            # get Ados dictionary
            N_he, he2idx, idx2he = Ados_dictionary(dims, tier)

            # start to construct the matrix L_he
            print("Start constructing process...")
            L_he = spzeros(ComplexF64, N_he * sup_dim, N_he * sup_dim)

            if liouville || progressBar
                print("\n")
            end
            flush(stdout)
            if progressBar
                prog = Progress(N_he; desc="Construct hierarchy matrix: ", PROGBAR_OPTIONS...)
            end
            for idx in 1:N_he
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
                L_he += pad_csc(- sum_ω * I_sup, N_he, N_he, idx, idx)

                state_neigh = copy(state)
                for k in 1:N_exp_term
                    n_k = state[k]
                    if n_k >= 1
                        state_neigh[k] = n_k - 1
                        idx_neigh = he2idx[state_neigh]
                        op = -1im * n_k * (η_list[k] * spreQ - conj(η_list[k]) * spostQ)
                        L_he += pad_csc(op, N_he, N_he, idx, idx_neigh)
                        state_neigh[k] = n_k
                    end
                    if n_exc <= tier - 1
                        state_neigh[k] = n_k + 1
                        idx_neigh = he2idx[state_neigh]
                        op = -1im * commQ
                        L_he += pad_csc(op, N_he, N_he, idx, idx_neigh)
                        state_neigh[k] = n_k
                    end
                end
                if progressBar
                    next!(prog)
                end
            end

            if liouville
                printstyled("Construct Liouvillian...", color=:green)
                flush(stdout)
                L_he += kron(sparse(I, N_he, N_he), liouvillian(Hsys, Jump_Ops, progressBar))
            end
            
            println("[DONE]\n")
            flush(stdout)
            return new(L_he, tier, Nsys, N_he, sup_dim)
        end
    end
    
    """
    # `evolution(M, ρ0, tlist; [solver, returnAdos, progressBar, SOLVEROptions...])`
    Solve the evolution (ODE problem) using HEOM model.

    ## Parameters
    - `M::AbstractHEOMMatrix` : the matrix given from HEOM model
    - `ρ0::Union{AbstractMatrix, AbstractOperator}` : initial state (density matrix)
    - `tlist::AbstractVector` : Denote the specific time points to save the solution at, during the solving process.
    - `solver` : solver in package `DifferentialEquations.jl`. Default to `Tsit5()`.
    - `returnAdos::Bool` : Decide to return the entire Ados vector or not. Default as `false`.
    - `progressBar::Bool` : Display progress bar during the process or not. Defaults to `true`.
    - `SOLVEROptions` : extra options for solver 

    ## Returns
    - `ρ_list` : The reduced density matrices in each time point.
    - `Ados_vec` : The Ados vector in each time point (only return when `returnAdos=true`).
    """
    function evolution(
            M::AbstractHEOMMatrix, 
            ρ0::Union{AbstractMatrix, AbstractOperator}, 
            tlist::AbstractVector;
            solver=Tsit5(),
            returnAdos::Bool=false,
            progressBar::Bool=true,
            SOLVEROptions...
        )
    
        (N1, N2) = size(ρ0)
        if (N1 != M.N_sys) || (N2 != M.N_sys)
            error("The size of initial state ρ0 is incorrect.")
        end

        # setup ρ_he and ρlist
        ρlist::Vector{SparseMatrixCSC{ComplexF64, Int64}} = []
        ρ_he::SparseVector{ComplexF64, Int64} = spzeros(M.N_he * M.sup_dim)
        if typeof(ρ0) <: AbstractMatrix
            push!(ρlist, ρ0)
            ρ_he[1:(M.sup_dim)] .= vec(ρ0)
        else
            push!(ρlist, ρ0.data)
            ρ_he[1:(M.sup_dim)] .= vec(ρ0.data)
        end
    
        if returnAdos
            Ados = [ρ_he]
        end
        
        # setup integrator
        dt_list = diff(tlist)
        integrator = init(
                ODEProblem(hierachy!, ρ_he, (tlist[1], tlist[end]), M.data),
                solver;
                SOLVEROptions...
            )
        
        # start solving ode
        print("Start solving hierachy equations of motions...")
        if progressBar
            print("\n")
            prog = Progress(length(tlist); start=1, desc="Progress : ", PROGBAR_OPTIONS...)
        end
        flush(stdout)
        for dt in dt_list
            step!(integrator, dt, true)
            
            # save the reduced density matrix and current ados vector (if \'returnAdos\' is specified)
            push!(ρlist, reshape(integrator.u[1:(M.sup_dim)], M.N_sys, M.N_sys))
            if returnAdos
                push!(Ados, integrator.u)
            end
        
            if progressBar
                next!(prog)
            end
        end
    
        println("[DONE]\n")
        flush(stdout)

        if returnAdos
            return ρlist, Ados
        else
            return ρlist
        end
    end

    # func. for solving evolution ODE
    function hierachy!(dρ, ρ, L, t)
        dρ .= L * ρ
    end
end
function A_minus(ω::Real, A::AbstractMatrix, iden::AbstractMatrix)
    return -1im * ω * iden + A
end

function A_plus(ω::Real, A::AbstractMatrix, iden::AbstractMatrix)
    return  1im * ω * iden + A
end

"""
# `DOS(M, ρ, ω_list, OP; [solver, progressBar, SOLVEROptions...])`
Calculate density of states.

## Parameters
- `M::AbstractHEOMMatrix` : the matrix given from HEOM model (the parity must be \":odd\".)
- `ρ::Union{AbstractMatrix, ADOs}` :  the system density matrix or the auxiliary density operators.
- `ω_list::AbstractVector` : the specific frequency points to solve.
- `OP::AbstractMatrix` : The system operator for the two-time correlation function in frequency domain.
- `solver` : solver in package `LinearSolve.jl`. Default to `UMFPACKFactorization()`.
- `progressBar::Bool` : Display progress bar during the process or not. Defaults to `true`.
- `SOLVEROptions` : extra options for solver 

## Returns
- `dos::AbstractVector` : density of state
"""
function DOS(
        M::AbstractHEOMMatrix, 
        ρ::T, 
        ω_list::AbstractVector, 
        OP::AbstractMatrix; 
        solver=UMFPACKFactorization(), 
        progressBar::Bool = true,
        SOLVEROptions...
    ) where T <: Union{AbstractMatrix, ADOs}

    # check parity
    if (M.parity != :odd)
        error("The parity of M must be \":odd\".")
    end    

    local b::AbstractVector

    # check ρ
    if T <: AbstractMatrix
        if size(ρ) == (M.dim, M.dim) 
            b = sparse(sparsevec(ρ))
        else
            error("The dimension of ρ should be equal to \"($(M.dim), $(M.dim))\".")
        end

    else # ρ::ADOs
        if (M.dim != ρ.dim)
            error("The system dimension between M and ρ are not consistent.")
        end
    
        if (M.Nb != ρ.Nb)
            error("The number of bosonic states between M and ρ are not consistent.")
        end
    
        if (M.Nf != ρ.Nf)
            error("The number of fermionic states between M and ρ are not consistent.")
        end

        b = ρ.data
    end

    # check dimension of OP
    if size(OP) != (M.dim, M.dim)
        error("The dimension of OP should be equal to \"($(M.dim), $(M.dim))\".")
    end

    Size, = size(M)
    I_total = sparse(I, Size, Size)
    I_heom  = sparse(I, M.N, M.N)

    # equal to : transpose(sparse(vec(system_identity_matrix)))
    I_dual_vec = transpose(sparsevec([1 + n * (M.dim + 1) for n in 0:(M.dim - 1)], ones(M.dim), M.sup_dim))

    # operators for calculating two-time correlation functions in frequency domain
    C_normal = kron(I_heom, spre(OP))
    C_dagger = kron(I_heom, spre(OP'))
    local b_minus::Vector{ComplexF64} = -1 * C_normal * b
    local b_plus ::Vector{ComplexF64} = -1 * C_dagger * b

    print("Start calculating density of states...")
    
    # solve for the first ω in ω_list to obtain the cache from the LinearSolve solution
    if progressBar
        print("\n")
        prog = Progress(length(ω_list); start=1, desc="Progress : ", PROGBAR_OPTIONS...)
    end
    ω = ω_list[1]
    sol_m = solve(LinearProblem(A_minus(ω, M.data, I_total), b_minus), solver, SOLVEROptions...)
    sol_p = solve(LinearProblem( A_plus(ω, M.data, I_total), b_plus ), solver, SOLVEROptions...)
    cache_m  = sol_m.cache
    cache_p  = sol_p.cache
    Cω_minus = C_dagger * sol_m.u
    Cω_plus  = C_normal * sol_p.u
    dos = [real(I_dual_vec * Cω_minus[1:(M.sup_dim)]) + real(I_dual_vec * Cω_plus[1:(M.sup_dim)])]
    if progressBar
        next!(prog)
    end 

    # take the cache to solve for the rest of the ω in ω_list
    for ω in ω_list[2:end]

        sol_m    = solve(set_A(cache_m, A_minus(ω, M.data, I_total)), solver, SOLVEROptions...)
        Cω_minus = C_dagger * sol_m.u
        
        sol_p    = solve(set_A(cache_p,  A_plus(ω, M.data, I_total)), solver, SOLVEROptions...)
        Cω_plus  = C_normal * sol_p.u

        # trace over the hilbert space of system (expectation value)
        push!(dos, real(I_dual_vec * Cω_minus[1:(M.sup_dim)]) + real(I_dual_vec * Cω_plus[1:(M.sup_dim)]))

        if progressBar
            next!(prog)
        end 
    end

    return dos
end
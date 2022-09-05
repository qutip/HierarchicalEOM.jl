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
- `M::AbstractHEOMMatrix` : the matrix given from HEOM model (the parity must be `:odd`.)
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

    Size, = size(M)
    I_total = sparse(I, Size, Size)
    I_heom  = sparse(I, M.N, M.N)

    local b::AbstractVector

    # check ρ
    if T <: AbstractMatrix
        if size(ρ) == (M.dim, M.dim) 
            v = sparsevec(ρ)
            b = sparsevec(v.nzind, v.nzval, Size)
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

    # equal to : transpose(sparse(vec(system_identity_matrix)))
    I_dual_vec = transpose(sparsevec([1 + n * (M.dim + 1) for n in 0:(M.dim - 1)], ones(M.dim), M.sup_dim))

    # operators for calculating two-time correlation functions in frequency domain
    C_normal = kron(I_heom, spre(OP))
    C_dagger = kron(I_heom, spre(OP'))
    local b_minus::Vector{ComplexF64} = -1 * C_normal * b
    local b_plus ::Vector{ComplexF64} = -1 * C_dagger * b

    print("Start calculating density of states...")
    Length = length(ω_list)
    dos    = Vector{Float64}(undef, length(Length))
    if progressBar
        print("\n")
        prog = Progress(Length; start=1, desc="Progress : ", PROGBAR_OPTIONS...)
    end
    @inbounds for i in 1:Length
        sol_m = solve(LinearProblem(A_minus(ω_list[i], M.data, I_total), b_minus), solver, SOLVEROptions...)
        sol_p = solve(LinearProblem( A_plus(ω_list[i], M.data, I_total), b_plus ), solver, SOLVEROptions...)
        Cω_minus = C_dagger * sol_m.u
        Cω_plus  = C_normal * sol_p.u

        # trace over the Hilbert space of system (expectation value)
        dos[i] = real(I_dual_vec * Cω_minus[1:(M.sup_dim)]) + real(I_dual_vec * Cω_plus[1:(M.sup_dim)])

        if progressBar
            next!(prog)
        end 
    end

    return dos
end
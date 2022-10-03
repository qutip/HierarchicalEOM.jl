"""
    DOS(M, ρ, ω_list, op; solver, progressBar, filename, SOLVEROptions...)
Calculate density of states.

# Parameters
- `M::AbstractHEOMMatrix` : the matrix given from HEOM model (the parity must be `:odd`.)
- `ρ` :  the system density matrix or the auxiliary density operators.
- `ω_list::AbstractVector` : the specific frequency points to solve.
- `op` : The system operator for the two-time correlation function in frequency domain.
- `solver` : solver in package `LinearSolve.jl`. Default to `UMFPACKFactorization()`.
- `progressBar::Bool` : Display progress bar during the process or not. Defaults to `true`.
- `filename::String` : If filename was specified, the value of dos for each ω will be saved into the file during the solving process.
- `SOLVEROptions` : extra options for solver 

For more details about solvers and extra options, please refer to [`LinearSolve.jl`](http://linearsolve.sciml.ai/stable/)

# Returns
- `dos::AbstractVector` : density of state
"""
function DOS(
        M::AbstractHEOMMatrix, 
        ρ, 
        ω_list::AbstractVector, 
        op; 
        solver=UMFPACKFactorization(), 
        progressBar::Bool = true,
        filename::String = "",
        SOLVEROptions...
    )

    SAVE::Bool = (filename != "")
    if SAVE && isfile(filename)
        error("FILE: $(filename) already exist.")
    end

    # check number of fermion states
    if M.Nf <= 0
        error("The number of fermionic states must be greater than zero, i.e., \"M.Nf > 0\".")

    # check parity
    elseif M.parity != :odd
        error("The parity of M must be \":odd\".")
    end

    Size, = size(M)
    I_total = sparse(I, Size, Size)
    I_heom  = sparse(I, M.N, M.N)

    local b::AbstractVector

    # check ρ
    if typeof(ρ) == ADOs  # ρ::ADOs
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
        
    elseif isValidMatrixType(ρ, M.dim)
        v = sparsevec(ρ)
        b = sparsevec(v.nzind, v.nzval, Size)
    else
        error("Invalid matrix \"ρ\".")
    end

    # check dimension of op
    if !isValidMatrixType(op, M.dim)
        error("Invalid matrix \"op\".")
    end

    # equal to : transpose(sparse(vec(system_identity_matrix)))
    I_dual_vec = transpose(sparsevec([1 + n * (M.dim + 1) for n in 0:(M.dim - 1)], ones(M.dim), M.sup_dim))

    # operators for calculating two-time correlation functions in frequency domain
    C_normal = kron(I_heom, spre(op))
    C_dagger = kron(I_heom, spre(op'))
    local b_minus::Vector{ComplexF64} = -1 * C_normal * b
    local b_plus ::Vector{ComplexF64} = -1 * C_dagger * b

    Length = length(ω_list)
    dos    = Vector{Float64}(undef, Length)
    print("Calculating density of states...")
    flush(stdout)
    if progressBar
        print("\n")
        prog = Progress(Length; start=1, desc="Progress : ", PROGBAR_OPTIONS...)
    end
    @inbounds for (i, ω) in enumerate(ω_list)
        if ω == 0
            sol_m = solve(LinearProblem(M.data, b_minus), solver, SOLVEROptions...)
            sol_p = solve(LinearProblem(M.data, b_plus),  solver, SOLVEROptions...)
        else
            Iω = 1im * ω * I_total
            sol_m = solve(LinearProblem(M.data - Iω, b_minus), solver, SOLVEROptions...)
            sol_p = solve(LinearProblem(M.data + Iω, b_plus),  solver, SOLVEROptions...)
        end
        Cω_minus = C_dagger * sol_m.u
        Cω_plus  = C_normal * sol_p.u
        
        # trace over the Hilbert space of system (expectation value)
        dos[i] = real(I_dual_vec * Cω_minus[1:(M.sup_dim)]) + real(I_dual_vec * Cω_plus[1:(M.sup_dim)])

        if SAVE
            open(filename, "a") do file
                write(file, "$(dos[i]),\n")
            end
        end

        if progressBar
            next!(prog)
        end 
    end
    println("[DONE]")

    return dos
end
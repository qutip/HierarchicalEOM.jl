"""
    PSD(M, ρ, op, ω_list; solver, verbose, filename, SOLVEROptions...)
Calculate power spectral density.

# Parameters
- `M::AbstractHEOMMatrix` : the matrix given from HEOM model (the parity must be `:none` or `:even`.)
- `ρ` :  the system density matrix or the auxiliary density operators.
- `op` : The system operator for the two-time correlation function in frequency domain.
- `ω_list::AbstractVector` : the specific frequency points to solve.
- `solver` : solver in package `LinearSolve.jl`. Default to `UMFPACKFactorization()`.
- `verbose::Bool` : To display verbose output and progress bar during the process or not. Defaults to `true`.
- `filename::String` : If filename was specified, the value of psd for each ω will be saved into the file during the solving process.
- `SOLVEROptions` : extra options for solver 

For more details about solvers and extra options, please refer to [`LinearSolve.jl`](http://linearsolve.sciml.ai/stable/)

# Returns
- `psd::AbstractVector` : power spectral density
"""
function PSD(
        M::AbstractHEOMMatrix, 
        ρ, 
        op, 
        ω_list::AbstractVector; 
        solver=UMFPACKFactorization(), 
        verbose::Bool = true,
        filename::String = "",
        SOLVEROptions...
    )

    SAVE::Bool = (filename != "")
    if SAVE && isfile(filename)
        error("FILE: $(filename) already exist.")
    end

    # check parity
    if (M.parity == :odd) || (typeof(M) == M_Fermion)
        error("The parity of M must be \":none\" (bosonic) or \":even\" (mixed) bath.")
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
        if (M.N != ρ.N)
            error("The number N between M and ρ are not consistent.")
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

    # operator for calculating two-time correlation functions in frequency domain
    C_normal = kron(I_heom, spre(op))
    C_dagger = kron(I_heom, spre(op'))
    local Cb::Vector{ComplexF64} = -1 * C_normal * b

    Length = length(ω_list)
    psd    = Vector{Float64}(undef, Length)
    
    if verbose
        print("Calculating power spectral density...\n")
        flush(stdout)
        prog = Progress(Length; start=1, desc="Progress : ", PROGBAR_OPTIONS...)
    end
    @inbounds for (i, ω) in enumerate(ω_list)
        if ω == 0
            sol = solve(LinearProblem(M.data, Cb), solver, SOLVEROptions...)
        else
            Iω = 1im * ω * I_total
            sol = solve(LinearProblem(M.data - Iω, Cb), solver, SOLVEROptions...)
        end
        Cω = C_dagger * sol.u

        # trace over the Hilbert space of system (expectation value)
        psd[i] = real(I_dual_vec * Cω[1:(M.sup_dim)])

        if SAVE
            open(filename, "a") do file
                write(file, "$(psd[i]),\n")
            end
        end

        if verbose
            next!(prog)
        end 
    end
    GC.gc()  # clean the garbage collector
    if verbose
        println("[DONE]")
    end
    
    return psd
end
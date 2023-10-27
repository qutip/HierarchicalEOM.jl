@doc raw"""
    spectrum(M, ρ, op, ω_list; solver, verbose, filename, SOLVEROptions...)
Calculate spectrum for the system.

# To calculate spectrum for bosonic systems (usually known as power spectrum):
```math
\pi S(\omega)=\textrm{Re}\left\{\int_0^\infty dt \langle A^\dagger(t) A(0)\rangle e^{-i\omega t}\right\},
```
remember to set the parameters: 
- `M::AbstractHEOMLSMatrix`: should be `EVEN` parity
- `op`: the operator ``A`` for bosonic system as shown above 

# To calculate spectrum for fermionic systems (usually known as density of states):
```math
    \pi A(\omega)=\textrm{Re}\left\{\int_0^\infty dt \left[\langle d(t) d^\dagger(0)\rangle^* + \langle d^\dagger(t) d(0)\rangle \right] e^{-i\omega t}\right\},
```
remember to set the parameters: 
- `M::AbstractHEOMLSMatrix`: should be `ODD` parity
- `op`: the (annihilation) operator ``d`` for fermionic system as shown above 

# Parameters
- `M::AbstractHEOMLSMatrix` : the matrix given from HEOM model.
- `ρ` :  the system density matrix or the auxiliary density operators.
- `op` : The annihilation operator acting on the system.
- `ω_list::AbstractVector` : the specific frequency points to solve.
- `solver` : solver in package `LinearSolve.jl`. Default to `UMFPACKFactorization()`.
- `verbose::Bool` : To display verbose output and progress bar during the process or not. Defaults to `true`.
- `filename::String` : If filename was specified, the value of spectrum for each ω will be saved into the file "filename.txt" during the solving process.
- `SOLVEROptions` : extra options for solver 

For more details about solvers and extra options, please refer to [`LinearSolve.jl`](http://linearsolve.sciml.ai/stable/)

# Returns
- `spec::AbstractVector` : the spectrum list corresponds to the specified `ω_list`
"""
function spectrum(
        M::AbstractHEOMLSMatrix, 
        ρ, 
        op, 
        ω_list::AbstractVector; 
        solver=UMFPACKFactorization(), 
        verbose::Bool = true,
        filename::String = "",
        SOLVEROptions...
    )

    Size = size(M, 1)

    # check ρ
    if typeof(ρ) == ADOs  # ρ::ADOs
        if (M.dim != ρ.dim)
            error("The system dimension between M and ρ are not consistent.")
        end
        if (M.N != ρ.N)
            error("The ADOs number \"N\" between M and ados are not consistent.")
        end
        if (typeof(ρ.parity) == OddParity)
            error("The parity of ρ or the ADOs must be `EVEN`.")
        end

        ados_vec = ρ.data
        
    else
        _ρ = HandleMatrixType(ρ, M.dim, "ρ (state)")
        v  = sparsevec(_ρ)
        ados_vec = sparsevec(v.nzind, v.nzval, Size)
    end

    # check dimension of op
    _op = HandleMatrixType(op, M.dim, "op (operator)")

    # check parity and calculate spectrum
    if (typeof(M.parity) == OddParity)
        return _density_of_states(M, ados_vec, _op, ω_list; 
            solver = solver, 
            verbose = verbose,
            filename = filename,
            SOLVEROptions...
        )
    else
        return _power_spectrum(M, ados_vec, _op, ω_list; 
            solver = solver, 
            verbose = verbose,
            filename = filename,
            SOLVEROptions...
        )
    end
end

@noinline function _power_spectrum(
        M::AbstractHEOMLSMatrix, 
        ados_vec::AbstractVector, 
        op,
        ω_list::AbstractVector; 
        solver=UMFPACKFactorization(), 
        verbose::Bool = true,
        filename::String = "",
        SOLVEROptions...
    )

    SAVE::Bool = (filename != "")
    if SAVE
        FILENAME = filename * ".txt"
        if isfile(FILENAME)
            error("FILE: $(FILENAME) already exist.")
        end
    end

    Size = size(M, 1)
    I_heom  = sparse(one(ComplexF64) * I, M.N, M.N)
    I_total = _HandleIdentityType(typeof(M.data), Size)

    # equal to : transpose(sparse(vec(system_identity_matrix)))
    I_dual_vec = transpose(sparsevec([1 + n * (M.dim + 1) for n in 0:(M.dim - 1)], ones(ComplexF64, M.dim), M.sup_dim))

    # operator for calculating two-time correlation functions in frequency domain
    a_normal = kron(I_heom, spre(op))
    a_dagger = kron(I_heom, spre(op'))
    X = _HandleVectorType(typeof(M.data), a_normal * ados_vec)

    Length = length(ω_list)
    Sω = Vector{Float64}(undef, Length)

    if verbose
        print("Calculating spectrum for bosonic systems...\n")
        flush(stdout)
        prog = Progress(Length; desc="Progress : ", PROGBAR_OPTIONS...)
    end
    Iω    = 1im * ω_list[1] * I_total
    cache = init(LinearProblem(M.data - Iω, X), solver, SOLVEROptions...)
    sol   = solve!(cache)
    @inbounds for (i, ω) in enumerate(ω_list)
        if i > 1            
            Iω  = 1im * ω * I_total
            cache.A = M.data - Iω
            sol = solve!(cache)
        end

        # trace over the Hilbert space of system (expectation value)
        Sω[i] = -1 * real(I_dual_vec * (a_dagger * _HandleVectorType(sol.u, false))[1:(M.sup_dim)])

        if SAVE
            open(FILENAME, "a") do file
                write(file, "$(Sω[i]),\n")
            end
        end

        if verbose
            next!(prog)
        end 
    end
    if verbose
        println("[DONE]")
    end

    return Sω
end

@noinline function _density_of_states(
        M::AbstractHEOMLSMatrix, 
        ados_vec::AbstractVector, 
        op,
        ω_list::AbstractVector; 
        solver=UMFPACKFactorization(), 
        verbose::Bool = true,
        filename::String = "",
        SOLVEROptions...
    )

    SAVE::Bool = (filename != "")
    if SAVE
        FILENAME = filename * ".txt"
        if isfile(FILENAME)
            error("FILE: $(FILENAME) already exist.")
        end
    end

    Size = size(M, 1)
    I_heom  = sparse(one(ComplexF64) * I, M.N, M.N)
    I_total = _HandleIdentityType(typeof(M.data), Size)

    # equal to : transpose(sparse(vec(system_identity_matrix)))
    I_dual_vec = transpose(sparsevec([1 + n * (M.dim + 1) for n in 0:(M.dim - 1)], ones(ComplexF64, M.dim), M.sup_dim))

    # operators for calculating two-time correlation functions in frequency domain
    d_normal = kron(I_heom, spre(op))
    d_dagger = kron(I_heom, spre(op'))
    X_m = _HandleVectorType(typeof(M.data), d_normal * ados_vec)
    X_p = _HandleVectorType(typeof(M.data), d_dagger * ados_vec)

    Length = length(ω_list)
    Aω = Vector{Float64}(undef, Length)

    if verbose
        print("Calculating spectrum for fermionic systems...\n")
        flush(stdout)
        prog = Progress(Length; desc="Progress : ", PROGBAR_OPTIONS...)
    end
    Iω = 1im * ω_list[1] * I_total
    cache_m = init(LinearProblem(M.data - Iω, X_m),  solver, SOLVEROptions...)
    cache_p = init(LinearProblem(M.data + Iω, X_p), solver, SOLVEROptions...)
    sol_m  = solve!(cache_m)
    sol_p  = solve!(cache_p)
    @inbounds for (i, ω) in enumerate(ω_list)
        if i > 1
            Iω = 1im * ω * I_total

            cache_m.A = M.data - Iω
            sol_m = solve!(cache_m)
            
            cache_p.A = M.data + Iω
            sol_p = solve!(cache_p)
        end
        Cω_m = d_dagger * _HandleVectorType(sol_m.u, false)
        Cω_p = d_normal * _HandleVectorType(sol_p.u, false)
        
        # trace over the Hilbert space of system (expectation value)
        Aω[i] = -1 * (
            real(I_dual_vec * Cω_p[1:(M.sup_dim)]) + 
            real(I_dual_vec * Cω_m[1:(M.sup_dim)])
        )

        if SAVE
            open(FILENAME, "a") do file
                write(file, "$(Aω[i]),\n")
            end
        end

        if verbose
            next!(prog)
        end 
    end
    if verbose
        println("[DONE]")
    end

    return Aω
end

function _HandleIdentityType(MatrixType::Type{TM}, S::Int) where TM <: SparseMatrixCSC
    ElType = eltype(MatrixType)
    return sparse(one(ElType) * I, S, S)
end
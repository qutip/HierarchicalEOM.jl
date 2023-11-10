@doc raw"""
    DensityOfStates(M, ρ, P_op, Q_op, ωlist; solver, verbose, filename, SOLVEROptions...)
Calculate density of states for the fermionic system.

```math
    \pi A(\omega)=\textrm{Re}\left\{\int_0^\infty dt \left[\langle d(t) d^\dagger(0)\rangle^* + \langle d^\dagger(t) d(0)\rangle \right] e^{-i\omega t}\right\},
```

# Parameters
- `M::AbstractHEOMLSMatrix` : the HEOMLS matrix which acts on `ODD`-parity operators.
- `ρ` :  the system density matrix or the auxiliary density operators.
- `op` : The annihilation operator (``d`` as shown above) acting on the fermionic system.
- `ωlist::AbstractVector` : the specific frequency points to solve.
- `solver` : solver in package `LinearSolve.jl`. Default to `UMFPACKFactorization()`.
- `verbose::Bool` : To display verbose output and progress bar during the process or not. Defaults to `true`.
- `filename::String` : If filename was specified, the value of spectrum for each ω will be saved into the file "filename.txt" during the solving process.
- `SOLVEROptions` : extra options for solver 

For more details about solvers and extra options, please refer to [`LinearSolve.jl`](http://linearsolve.sciml.ai/stable/)

# Returns
- `dos::AbstractVector` : the list of density of states corresponds to the specified `ωlist`
"""
@noinline function DensityOfStates(
        M::AbstractHEOMLSMatrix, 
        ρ,
        op,
        ωlist::AbstractVector; 
        solver=UMFPACKFactorization(), 
        verbose::Bool = true,
        filename::String = "",
        SOLVEROptions...
    )

    Size = size(M, 1)

    # check M
    if (typeof(M.parity) != OddParity)
        error("The HEOMLS matrix M must be acting on `ODD`-parity operators.")
    end

    # check ρ
    if typeof(ρ) == ADOs  # ρ::ADOs
        
        _check_sys_dim_and_ADOs_num(M, ρ)

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

    SAVE::Bool = (filename != "")
    if SAVE
        FILENAME = filename * ".txt"
        if isfile(FILENAME)
            error("FILE: $(FILENAME) already exist.")
        end
    end

    I_heom  = sparse(one(ComplexF64) * I, M.N, M.N)
    I_total = _HandleIdentityType(typeof(M.data), Size)

    # equal to : transpose(sparse(vec(system_identity_matrix)))
    I_dual_vec = transpose(sparsevec([1 + n * (M.dim + 1) for n in 0:(M.dim - 1)], ones(ComplexF64, M.dim), M.sup_dim))

    # operators for calculating two-time correlation functions in frequency domain
    d_normal = kron(I_heom, spre(_op))
    d_dagger = kron(I_heom, spre(_op'))
    X_m = _HandleVectorType(typeof(M.data), d_normal * ados_vec)
    X_p = _HandleVectorType(typeof(M.data), d_dagger * ados_vec)

    ElType = eltype(M)
    ωList  = _HandleFloatType(ElType, ωlist)
    Length = length(ωList)
    Aω = Vector{Float64}(undef, Length)

    if verbose
        print("Calculating density of states...\n")
        flush(stdout)
        prog = Progress(Length; desc="Progress : ", PROGBAR_OPTIONS...)
    end
    i  = convert(ElType, 1im)
    Iω = i * ωList[1] * I_total
    cache_m = init(LinearProblem(M.data - Iω, X_m),  solver, SOLVEROptions...)
    cache_p = init(LinearProblem(M.data + Iω, X_p), solver, SOLVEROptions...)
    sol_m  = solve!(cache_m)
    sol_p  = solve!(cache_p)
    @inbounds for (j, ω) in enumerate(ωList)
        if j > 1
            Iω = i * ω * I_total

            cache_m.A = M.data - Iω
            sol_m = solve!(cache_m)
            
            cache_p.A = M.data + Iω
            sol_p = solve!(cache_p)
        end
        Cω_m = d_dagger * _HandleVectorType(sol_m.u, false)
        Cω_p = d_normal * _HandleVectorType(sol_p.u, false)
        
        # trace over the Hilbert space of system (expectation value)
        Aω[j] = -1 * (
            real(I_dual_vec * Cω_p[1:(M.sup_dim)]) + 
            real(I_dual_vec * Cω_m[1:(M.sup_dim)])
        )

        if SAVE
            open(FILENAME, "a") do file
                write(file, "$(Aω[j]),\n")
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
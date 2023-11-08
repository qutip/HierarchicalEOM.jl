@doc raw"""
    spectrum(M, ρ, op, ωlist; solver, verbose, filename, SOLVEROptions...)
!!! warning "Warning"
    This function has been deprecated start from `HierarchicalEOM v1.1`, use `PowerSpectrum` or `DensityOfStates` instead.
"""
function spectrum(
        M::AbstractHEOMLSMatrix, 
        ρ, 
        op, 
        ωlist::AbstractVector; 
        solver=UMFPACKFactorization(), 
        verbose::Bool = true,
        filename::String = "",
        SOLVEROptions...
    )
    error("This function has been deprecated start from \`HierarchicalEOM v1.1\`, use \`PowerSpectrum\` or \`DensityOfStates\` instead.")
end

@doc raw"""
    PowerSpectrum(M, ρ, Q_op, ωlist, reverse; solver, verbose, filename, SOLVEROptions...)
Calculate spectrum for the system where `P_op` will be automatically set as the adjoint of `Q_op`.

This function is equivalent to:
`PowerSpectrum(M, ρ, Q_op', Q_op, ωlist, reverse; solver, verbose, filename, SOLVEROptions...)`
"""
function PowerSpectrum(
        M::AbstractHEOMLSMatrix, 
        ρ, 
        Q_op, 
        ωlist::AbstractVector,
        reverse::Bool = false; 
        solver=UMFPACKFactorization(), 
        verbose::Bool = true,
        filename::String = "",
        SOLVEROptions...
    )
    return PowerSpectrum(M, ρ, Q_op', Q_op, ωlist, reverse;
        solver = solver, 
        verbose = verbose,
        filename = filename,
        SOLVEROptions...
    )
end

@doc raw"""
    PowerSpectrum(M, ρ, P_op, Q_op, ωlist, reverse; solver, verbose, filename, SOLVEROptions...)
Calculate power spectrum for the system.

```math
\pi S(\omega)=\textrm{Re}\left\{\int_0^\infty dt \langle P(t) Q(0)\rangle e^{-i\omega t}\right\},
```

# To calculate spectrum when input operator `Q_op` has `EVEN`-parity:
remember to set the parameters: 
- `M::AbstractHEOMLSMatrix`: should be `EVEN` parity

# To calculate spectrum when input operator `Q_op` has `ODD`-parity:
remember to set the parameters: 
- `M::AbstractHEOMLSMatrix`: should be `ODD` parity

# Parameters
- `M::AbstractHEOMLSMatrix` : the HEOMLS matrix.
- `ρ` :  the system density matrix or the auxiliary density operators.
- `P_op`: the operator ``P`` acting on the system.
- `Q_op`: the operator ``Q`` acting on the system.
- `ωlist::AbstractVector` : the specific frequency points to solve.
- `reverse::Bool` : If `true`, calculate ``\langle P(-t)Q(0) \rangle = \langle P(0)Q(t) \rangle = \langle P(t)Q(0) \rangle^*`` instead of ``\langle P(t) Q(0) \rangle``. Default to `false`.
- `solver` : solver in package `LinearSolve.jl`. Default to `UMFPACKFactorization()`.
- `verbose::Bool` : To display verbose output and progress bar during the process or not. Defaults to `true`.
- `filename::String` : If filename was specified, the value of spectrum for each ω will be saved into the file "filename.txt" during the solving process.
- `SOLVEROptions` : extra options for solver 

For more details about solvers and extra options, please refer to [`LinearSolve.jl`](http://linearsolve.sciml.ai/stable/)

# Returns
- `spec::AbstractVector` : the spectrum list corresponds to the specified `ωlist`
"""
function PowerSpectrum(
        M::AbstractHEOMLSMatrix, 
        ρ, 
        P_op,
        Q_op, 
        ωlist::AbstractVector,
        reverse::Bool = false; 
        solver = UMFPACKFactorization(), 
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

    # check dimension of P_op and Q_op
    _P = HandleMatrixType(P_op, M.dim, "P_op (operator)")
    _Q = HandleMatrixType(Q_op, M.dim, "Q_op (operator)")

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

    # operator for calculating two-time correlation functions in frequency domain
    P_sup = kron(I_heom, spre(_P))
    Q_sup = kron(I_heom, spre(_Q))
    X = _HandleVectorType(typeof(M.data), Q_sup * ados_vec)

    ElType = eltype(M)
    ωList  = _HandleFloatType(ElType, ωlist)
    Length = length(ωList)
    Sω = Vector{Float64}(undef, Length)

    if verbose
        print("Calculating spectrum (with EVEN-parity operators)...\n")
        flush(stdout)
        prog = Progress(Length; desc="Progress : ", PROGBAR_OPTIONS...)
    end
    if reverse
        i = convert(ElType,  1im)
    else
        i = convert(ElType, -1im)
    end
    Iω    = i * ωList[1] * I_total
    cache = init(LinearProblem(M.data + Iω, X), solver, SOLVEROptions...)
    sol   = solve!(cache)
    @inbounds for (j, ω) in enumerate(ωList)
        if j > 1            
            Iω  = i * ω * I_total
            cache.A = M.data + Iω
            sol = solve!(cache)
        end

        # trace over the Hilbert space of system (expectation value)
        Sω[j] = -1 * real(I_dual_vec * (P_sup * _HandleVectorType(sol.u, false))[1:(M.sup_dim)])

        if SAVE
            open(FILENAME, "a") do file
                write(file, "$(Sω[j]),\n")
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

function _HandleIdentityType(MatrixType::Type{TM}, S::Int) where TM <: SparseMatrixCSC
    ElType = eltype(MatrixType)
    return sparse(one(ElType) * I, S, S)
end
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

function _HandleIdentityType(MatrixType::Type{TM}, S::Int) where TM <: SparseMatrixCSC
    ElType = eltype(MatrixType)
    return sparse(one(ElType) * I, S, S)
end

@doc raw"""
    PowerSpectrum(M, ρ, Q_op, ωlist, reverse; solver, verbose, filename, SOLVEROptions...)
Calculate power spectrum for the system where `P_op` will be automatically set as the adjoint of `Q_op`.

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
        print("Calculating power spectrum...\n")
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
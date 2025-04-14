export DensityOfStates

@doc raw"""
    DensityOfStates(M, ρ, d_op, ωlist; solver, verbose, filename, SOLVEROptions...)
Calculate density of states for the fermionic system in frequency domain.

```math
    \pi A(\omega)=\textrm{Re}\left\{\int_0^\infty dt \left[\langle d(t) d^\dagger(0)\rangle^* + \langle d^\dagger(t) d(0)\rangle \right] e^{-i\omega t}\right\},
```

# Parameters
- `M::AbstractHEOMLSMatrix` : the HEOMLS matrix which acts on `ODD`-parity operators.
- `ρ::Union{QuantumObject,ADOs}` :  the system density matrix or the auxiliary density operators.
- `d_op::QuantumObject` : The annihilation operator (``d`` as shown above) acting on the fermionic system.
- `ωlist::AbstractVector` : the specific frequency points to solve.
- `solver::SciMLLinearSolveAlgorithm` : solver in package `LinearSolve.jl`. Default to `KrylovJL_BICGSTAB(rtol=1e-12, atol=1e-14)`.
- `verbose::Bool` : To display verbose output and progress bar during the process or not. Defaults to `true`.
- `filename::String` : If filename was specified, the value of spectrum for each ω will be saved into the file "filename.txt" during the solving process.
- `SOLVEROptions` : extra options for solver 

# Notes
- For more details about `solver` and `SOLVEROptions`, please refer to [`LinearSolve.jl`](http://linearsolve.sciml.ai/stable/)

# Returns
- `dos::AbstractVector` : the list of density of states corresponds to the specified `ωlist`
"""
@noinline function DensityOfStates(
    M::AbstractHEOMLSMatrix{<:MatrixOperator},
    ρ::Union{QuantumObject,ADOs},
    d_op::QuantumObject,
    ωlist::AbstractVector;
    solver::SciMLLinearSolveAlgorithm = KrylovJL_BICGSTAB(rtol = 1e-12, atol = 1e-14),
    verbose::Bool = true,
    filename::String = "",
    SOLVEROptions...,
)

    # check M
    if M.parity == EVEN
        error("The HEOMLS matrix M must be acting on `ODD`-parity operators.")
    end

    # Handle ρ
    if typeof(ρ) == ADOs  # ρ::ADOs
        ados = ρ
        if ados.parity != EVEN
            error("The parity of ρ must be `EVEN`.")
        end
    else
        ados = ADOs(ρ, M.N)
    end
    _check_sys_dim_and_ADOs_num(M, ados)

    # Handle d_op
    _tr = _Tr(M)
    Id_sys = I(prod(d_op.dimensions))
    Id_HEOM = I(M.N)
    d_normal = HEOMSuperOp(spre(d_op, Id_sys), ODD, M; Id_cache = Id_HEOM)
    d_dagger = HEOMSuperOp(spre(d_op', Id_sys), ODD, M; Id_cache = Id_HEOM)
    b_m = _HandleVectorType(M, (d_normal * ados).data)
    b_p = _HandleVectorType(M, (d_dagger * ados).data)
    _tr_d_normal = _HandleTraceVectorType(M, adjoint(d_normal.data) * _tr) # another adjoint will be applied in dot function later
    _tr_d_dagger = _HandleTraceVectorType(M, adjoint(d_dagger.data) * _tr) # another adjoint will be applied in dot function later

    SAVE::Bool = (filename != "")
    if SAVE
        FILENAME = filename * ".txt"
        isfile(FILENAME) && error("FILE: $(FILENAME) already exist.")
    end

    ElType = eltype(M)
    ωList = convert(Vector{_FType(M)}, ωlist) # Convert it to support GPUs and avoid type instabilities
    Length = length(ωList)
    Aω = Vector{Float64}(undef, Length)

    if verbose
        print("Calculating density of states in frequency domain...\n")
        flush(stdout)
    end
    prog = ProgressBar(Length; enable = verbose)
    i = convert(ElType, 1im)
    I_total = I(size(M, 1))
    cache_m = cache_p = nothing
    for ω in ωList
        Iω = i * ω * I_total

        if prog.counter[] == 0
            cache_m = init(LinearProblem(M.data.A - Iω, b_m), solver, SOLVEROptions...)
            sol_m = solve!(cache_m)

            cache_p = init(LinearProblem(M.data.A + Iω, b_p), solver, SOLVEROptions...)
            sol_p = solve!(cache_p)
        else
            cache_m.A = M.data.A - Iω
            sol_m = solve!(cache_m)

            cache_p.A = M.data.A + Iω
            sol_p = solve!(cache_p)
        end

        # trace over the Hilbert space of system (expectation value)
        val = -1 * real(dot(_tr_d_normal, sol_p.u) + dot(_tr_d_dagger, sol_m.u))
        Aω[prog.counter[]+1] = val

        if SAVE
            open(FILENAME, "a") do file
                return write(file, "$(val),\n")
            end
        end
        next!(prog)
    end
    if verbose
        println("[DONE]")
    end

    return Aω
end

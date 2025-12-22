export DensityOfStates

@doc raw"""
    DensityOfStates(M, ρ, d_op, ωlist; alg, progress_bar, filename, kwargs...)
Calculate density of states for the fermionic system in frequency domain.

```math
    \pi A(\omega)=\textrm{Re}\left\{\int_0^\infty dt \left[\langle d(t) d^\dagger(0)\rangle^* + \langle d^\dagger(t) d(0)\rangle \right] e^{-i\omega t}\right\},
```

# Parameters
- `M::AbstractHEOMLSMatrix` : the HEOMLS matrix which acts on `ODD`-parity operators.
- `ρ::Union{QuantumObject,ADOs}` :  the system density matrix or the auxiliary density operators.
- `d_op::QuantumObject` : The annihilation operator (``d`` as shown above) acting on the fermionic system.
- `ωlist::AbstractVector` : the specific frequency points to solve.
- `alg::SciMLLinearSolveAlgorithm` : The solving algorithm in package `LinearSolve.jl`. Default to `KrylovJL_GMRES(rtol=1e-12, atol=1e-14)`.
- `progress_bar::Union{Val,Bool}`: Whether to show the progress bar. Defaults to `Val(true)`.
- `filename::String` : If filename was specified, the value of spectrum for each ω will be saved into the file "filename.txt" during the solving process.
- `kwargs` : The keyword arguments for `LinearProblem`.

# Notes
- For more details about `alg`, `kwargs`, and `LinearProblem`, please refer to [`LinearSolve.jl`](http://linearsolve.sciml.ai/stable/)

# Returns
- `dos::AbstractVector` : the list of density of states corresponds to the specified `ωlist`
"""
@noinline function DensityOfStates(
    M::AbstractHEOMLSMatrix,
    ρ::Union{QuantumObject,ADOs},
    d_op::QuantumObject,
    ωlist::AbstractVector;
    alg::SciMLLinearSolveAlgorithm = KrylovJL_GMRES(rtol = 1e-12, atol = 1e-14),
    progress_bar::Union{Val,Bool} = Val(true),
    filename::String = "",
    kwargs...,
)
    isconstant(M) || throw(ArgumentError("The HEOMLS matrix M must be time-independent to calculate DensityOfStates."))
    haskey(kwargs, :solver) &&
        error("The keyword argument `solver` for DensityOfStates is deprecated, use `alg` instead.")
    haskey(kwargs, :verbose) &&
        error("The keyword argument `verbose` for DensityOfStates is deprecated, use `progress_bar` instead.")

    # check M
    (M.parity == EVEN) && error("The HEOMLS matrix M must be acting on `ODD`-parity operators.")

    # Handle ρ
    if typeof(ρ) == ADOs  # ρ::ADOs
        ados = ρ
        (ados.parity != EVEN) && error("The parity of ρ must be `EVEN`.")
    else
        ados = ADOs(ρ, M.N)
    end
    _check_sys_dim_and_ADOs_num(M, ados)

    # Handle d_op
    _tr = _Tr(M)
    d_normal = HEOMSuperOp(spre(d_op), ODD, M)
    d_dagger = HEOMSuperOp(spre(d_op'), ODD, M)
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
    ωList = convert(Vector{_float_type(M)}, ωlist) # Convert it to support GPUs and avoid type instabilities
    Length = length(ωList)
    Aω = Vector{Float64}(undef, Length)

    progr = Progress(
        Length;
        enabled = getVal(progress_bar),
        desc = "[DensityOfStates] ",
        QuantumToolbox.settings.ProgressMeterKWARGS...,
    )
    i = convert(ElType, 1im)
    I_total = Eye(size(M, 1))
    Iω1 = i * ωList[1] * I_total
    A0 = needs_concrete_A(alg) ? concretize(M.data) : cache_operator(M.data, b_m)
    cache_m = init(LinearProblem(A0 - Iω1, b_m), alg, kwargs...)
    cache_p = init(LinearProblem(A0 + Iω1, b_p), alg, kwargs...)
    for (idx, ω) in enumerate(ωList)
        if idx > 1
            Iω = i * ω * I_total
            cache_m.A = A0 - Iω
            cache_p.A = A0 + Iω
        end
        sol_m = solve!(cache_m)
        sol_p = solve!(cache_p)

        # trace over the Hilbert space of system (expectation value)
        val = -1 * real(dot(_tr_d_normal, sol_p.u) + dot(_tr_d_dagger, sol_m.u))
        Aω[idx] = val

        if SAVE
            open(FILENAME, "a") do file
                return write(file, "$(val),\n")
            end
        end
        next!(progr)
    end
    return Aω
end

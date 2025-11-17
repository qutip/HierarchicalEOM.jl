export PowerSpectrum

@doc raw"""
    PowerSpectrum(M, ρ, Q_op, ωlist, reverse; alg, progress_bar, filename, kwargs...)
Calculate power spectrum for the system in frequency domain where `P_op` will be automatically set as the adjoint of `Q_op`.

This function is equivalent to:
`PowerSpectrum(M, ρ, Q_op', Q_op, ωlist, reverse; alg, progress_bar, filename, kwargs...)`
"""
PowerSpectrum(
    M::AbstractHEOMLSMatrix{<:MatrixOperator},
    ρ::Union{QuantumObject,ADOs},
    Q_op::QuantumObject,
    ωlist::AbstractVector,
    reverse::Bool = false;
    alg::SciMLLinearSolveAlgorithm = KrylovJL_GMRES(rtol = 1e-12, atol = 1e-14),
    progress_bar::Union{Val,Bool} = Val(true),
    filename::String = "",
    kwargs...,
) = PowerSpectrum(
    M,
    ρ,
    Q_op',
    Q_op,
    ωlist,
    reverse;
    alg = alg,
    progress_bar = progress_bar,
    filename = filename,
    kwargs...,
)

@doc raw"""
    PowerSpectrum(M, ρ, P_op, Q_op, ωlist, reverse; alg, progress_bar, filename, kwargs...)
Calculate power spectrum for the system in frequency domain.

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
- `ρ::Union{QuantumObject,ADOs}` :  the system density matrix or the auxiliary density operators.
- `P_op::Union{QuantumObject,HEOMSuperOp}`: the system operator (or `HEOMSuperOp`) ``P`` acting on the system.
- `Q_op::Union{QuantumObject,HEOMSuperOp}`: the system operator (or `HEOMSuperOp`) ``Q`` acting on the system.
- `ωlist::AbstractVector` : the specific frequency points to solve.
- `reverse::Bool` : If `true`, calculate ``\langle P(-t)Q(0) \rangle = \langle P(0)Q(t) \rangle = \langle P(t)Q(0) \rangle^*`` instead of ``\langle P(t) Q(0) \rangle``. Default to `false`.
- `alg::SciMLLinearSolveAlgorithm` : The solving algorithm in package `LinearSolve.jl`. Default to `KrylovJL_GMRES(rtol=1e-12, atol=1e-14)`.
- `progress_bar::Union{Val,Bool}`: Whether to show the progress bar. Defaults to `Val(true)`.
- `filename::String` : If filename was specified, the value of spectrum for each ω will be saved into the file "filename.txt" during the solving process.
- `kwargs` : The keyword arguments for `LinearProblem`.

# Notes
- For more details about `alg`, `kwargs`, and `LinearProblem`, please refer to [`LinearSolve.jl`](http://linearsolve.sciml.ai/stable/)

# Returns
- `spec::AbstractVector` : the spectrum list corresponds to the specified `ωlist`
"""
@noinline function PowerSpectrum(
    M::AbstractHEOMLSMatrix{<:MatrixOperator},
    ρ::Union{QuantumObject,ADOs},
    P_op,
    Q_op,
    ωlist::AbstractVector,
    reverse::Bool = false;
    alg::SciMLLinearSolveAlgorithm = KrylovJL_GMRES(rtol = 1e-12, atol = 1e-14),
    progress_bar::Union{Val,Bool} = Val(true),
    filename::String = "",
    kwargs...,
)
    haskey(kwargs, :solver) &&
        error("The keyword argument `solver` for PowerSpectrum is deprecated, use `alg` instead.")
    haskey(kwargs, :verbose) &&
        error("The keyword argument `verbose` for PowerSpectrum is deprecated, use `progress_bar` instead.")

    # Handle ρ
    if ρ isa ADOs
        ados = ρ
    else
        ados = ADOs(ρ, M.N)
    end
    _check_sys_dim_and_ADOs_num(M, ados)

    # Handle P_op
    if P_op isa HEOMSuperOp
        _check_sys_dim_and_ADOs_num(M, P_op)
        _P = P_op
    else
        _P = HEOMSuperOp(spre(P_op), EVEN, M)
    end
    _tr_P = _HandleTraceVectorType(M, adjoint(_P.data) * _Tr(M)) # another adjoint will be applied in dot function later

    # Handle Q_op
    if Q_op isa HEOMSuperOp
        _check_sys_dim_and_ADOs_num(M, Q_op)
        _Q_ados = Q_op * ados
        _check_parity(M, _Q_ados)
    else
        if M.parity == EVEN
            new_parity = ados.parity
        else
            new_parity = !ados.parity
        end
        _Q_ados = HEOMSuperOp(spre(Q_op), new_parity, M) * ados
    end
    b = _HandleVectorType(M, _Q_ados.data)

    SAVE::Bool = (filename != "")
    if SAVE
        FILENAME = filename * ".txt"
        isfile(FILENAME) && error("FILE: $(FILENAME) already exist.")
    end

    ElType = eltype(M)
    ωList = convert(Vector{_float_type(M)}, ωlist) # Convert it to support GPUs and avoid type instabilities
    Length = length(ωList)
    Sω = Vector{Float64}(undef, Length)

    progr = Progress(
        Length;
        enabled = getVal(progress_bar),
        desc = "[PowerSpectrum] ",
        QuantumToolbox.settings.ProgressMeterKWARGS...,
    )
    i = reverse ? convert(ElType, 1im) : i = convert(ElType, -1im)
    I_total = Eye(size(M, 1))
    cache = init(LinearProblem(M.data.A + i * ωList[1] * I_total, b), alg, kwargs...)
    for (idx, ω) in enumerate(ωList)
        if idx > 1
            cache.A = M.data.A + i * ω * I_total
        end
        sol = solve!(cache)

        # trace over the Hilbert space of system (expectation value)
        val = -1 * real(dot(_tr_P, sol.u))
        Sω[idx] = val

        if SAVE
            open(FILENAME, "a") do file
                return write(file, "$(val),\n")
            end
        end
        next!(progr)
    end
    return Sω
end

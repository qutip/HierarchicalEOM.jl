@doc raw"""
    correlation_3op_2t(
        M::AbstractHEOMLSMatrix,
        state::Union{QuantumObject,ADOs},
        tlist::AbstractVector,
        τlist::AbstractVector,
        A::QuantumObject{Operator},
        B::QuantumObject{Operator},
        C::QuantumObject{Operator};
        alg::AbstractODEAlgorithm = DP5(),
        max_t_plus_tau::Union{Nothing,Real} = nothing,
        ensemblealg::EnsembleAlgorithm = EnsembleThreads(),
        progress_bar::Union{Val,Bool} = Val(true),
        kwargs...,
    )

Returns the two-times correlation function of three operators ``\hat{A}``, ``\hat{B}`` and ``\hat{C}``: ``\left\langle \hat{A}(t) \hat{B}(t + \tau) \hat{C}(t) \right\rangle`` for a given `state` using HEOM approach.

# Notes

This function is implemented based on [`HEOMsolve_map`](@ref) by the following adjustments:

- If `max_t_plus_tau` is specified, skip computation where ``t + \tau >`` `max_t_plus_tau`. Skipped entries are filled with ``0``. Default `nothing` means compute all entries (equivalent to `Inf`).
- The ``\tau``-evolution starting from every `t` in `tlist` is solved together as a single [`HEOMsolve_map`](@ref), and `ensemblealg` (default to `EnsembleThreads()`) controls how that ensemble is parallelized.
"""
function QuantumToolbox.correlation_3op_2t(
        M::AbstractHEOMLSMatrix,
        state::T_state,
        tlist::AbstractVector,
        τlist::AbstractVector,
        A::QuantumObject{Operator},
        B::QuantumObject{Operator},
        C::QuantumObject{Operator};
        alg::AbstractODEAlgorithm = DP5(),
        max_t_plus_tau::Union{Nothing, Real} = nothing,
        ensemblealg::EnsembleAlgorithm = EnsembleThreads(),
        progress_bar::Union{Val, Bool} = Val(true),
        kwargs...,
    ) where {T_state <: Union{QuantumObject, ADOs}}

    # check tlist and τlist
    QuantumToolbox._check_correlation_time_list(tlist)
    QuantumToolbox._check_correlation_time_list(τlist)
    max_time = something(max_t_plus_tau, Inf)
    max_time >= 0 || throw(ArgumentError("`max_t_plus_tau` must be `nothing` or a non-negative number."))

    # number of τ values with t + τ ≤ max_time for each t (same as QuTiP's searchsorted(..., side='right'))
    # tlist is sorted, so the t values without any τ left (t > max_time) are all at the end and are skipped
    nτ_list = [searchsortedlast(τlist, max_time - t) for t in tlist]
    tlist_kept = tlist[nτ_list .> 0]

    kwargs2 = merge((saveat = collect(tlist),), (; kwargs...))

    # note that doesn't need to do the first HEOMsolve if tlist_kept == [0] (for 2op_1t and 3op_1t cases)
    ados_t_list = (tlist_kept == [0]) ? [state] : HEOMsolve(M, state, tlist_kept; alg, progress_bar, kwargs2...).ados

    # every t is one trajectory of a single HEOMsolve_map ensemble, evolving C ρ(t) A over t + τ and recording ⟨B⟩
    C_A = HEOMSuperOp(sprepost(C, A), M.parity, M)
    probs = [
        HEOMsolveProblem(M, C_A * ados_t_list[i], τlist[1:nτ_list[i]] .+ t; e_ops = [B], progress_bar = Val(false), kwargs...)
            for (i, t) in enumerate(tlist_kept)
    ]

    # each problem already has its own e_ops callback, so only the HEOMLS matrix (its coefficients are updated in place) is copied
    _prob_func(_, ctx) = (prob = probs[ctx.sim_id].prob; remake(prob, f = _copy_for_trajectory(prob.f.f)))
    sols = HEOMsolve_map(first(probs), probs, alg, ensemblealg; prob_func = _prob_func, safetycopy = false, progress_bar)

    # make the output correlation Matrix align with QuTiP
    # 1st dimension corresponds to tlist
    # 2nd dimension corresponds to τlist (the entries with t + τ > max_t_plus_tau are left as 0)
    corr = zeros(eltype(first(sols).expect), length(tlist), length(τlist))
    for (i, sol) in enumerate(sols)
        corr[i, 1:nτ_list[i]] .= @view sol.expect[1, :]
    end
    return corr
end

@doc raw"""
    correlation_3op_1t(
        M::AbstractHEOMLSMatrix,
        state::Union{QuantumObject,ADOs},
        τlist::AbstractVector,
        A::QuantumObject{Operator},
        B::QuantumObject{Operator},
        C::QuantumObject{Operator};
        alg::AbstractODEAlgorithm = DP5(),
        progress_bar::Union{Val,Bool} = Val(true),
        kwargs...,
    )

Returns the two-time correlation function (with only one time coordinate ``\tau``) of three operators ``\hat{A}``, ``\hat{B}`` and ``\hat{C}``: ``\left\langle \hat{A}(0) \hat{B}(\tau) \hat{C}(0) \right\rangle`` for a given `state` using HEOM approach.
"""
function QuantumToolbox.correlation_3op_1t(
        M::AbstractHEOMLSMatrix,
        state::T_state,
        τlist::AbstractVector,
        A::QuantumObject{Operator},
        B::QuantumObject{Operator},
        C::QuantumObject{Operator};
        alg::AbstractODEAlgorithm = DP5(),
        progress_bar::Union{Val, Bool} = Val(true),
        kwargs...,
    ) where {T_state <: Union{QuantumObject, ADOs}}
    # there is only a single t = 0, so the underlying HEOMsolve_map ensemble has just one trajectory
    corr = correlation_3op_2t(M, state, [0], τlist, A, B, C; alg, ensemblealg = EnsembleSerial(), progress_bar, kwargs...)
    return corr[1, :]
end

@doc raw"""
    correlation_2op_2t(
        M::AbstractHEOMLSMatrix,
        state::Union{QuantumObject,ADOs},
        tlist::AbstractVector,
        τlist::AbstractVector,
        A::QuantumObject{Operator},
        B::QuantumObject{Operator};
        alg::AbstractODEAlgorithm = DP5(),
        reverse::Bool = false,
        max_t_plus_tau::Union{Nothing,Real} = nothing,
        ensemblealg::EnsembleAlgorithm = EnsembleThreads(),
        progress_bar::Union{Val,Bool} = Val(true),
        kwargs...,
    )

Returns the two-times correlation function of two operators ``\hat{A}`` and ``\hat{B}`` : ``\left\langle \hat{A}(t + \tau) \hat{B}(t) \right\rangle`` for a given `state` using HEOM approach.

# Notes

- When `reverse=true`, the correlation function is calculated as ``\left\langle \hat{A}(t) \hat{B}(t + \tau) \right\rangle``.
- If `max_t_plus_tau` is specified, skip computation where ``t + \tau >`` `max_t_plus_tau`. Skipped entries are filled with ``0``. Default `nothing` means compute all entries (equivalent to `Inf`).
- The ``\tau``-evolution starting from every `t` in `tlist` is solved together as a single [`HEOMsolve_map`](@ref), and `ensemblealg` (default to `EnsembleThreads()`) controls how that ensemble is parallelized.
"""
function QuantumToolbox.correlation_2op_2t(
        M::AbstractHEOMLSMatrix,
        state::T_state,
        tlist::AbstractVector,
        τlist::AbstractVector,
        A::QuantumObject{Operator},
        B::QuantumObject{Operator};
        reverse::Bool = false,
        max_t_plus_tau::Union{Nothing, Real} = nothing,
        alg::AbstractODEAlgorithm = DP5(),
        ensemblealg::EnsembleAlgorithm = EnsembleThreads(),
        progress_bar::Union{Val, Bool} = Val(true),
        kwargs...,
    ) where {T_state <: Union{QuantumObject, ADOs}}
    C = qeye_like(A)

    if reverse
        corr = correlation_3op_2t(M, state, tlist, τlist, A, B, C; alg, max_t_plus_tau, ensemblealg, progress_bar, kwargs...)
    else
        corr = correlation_3op_2t(M, state, tlist, τlist, C, A, B; alg, max_t_plus_tau, ensemblealg, progress_bar, kwargs...)
    end

    return corr
end

@doc raw"""
    correlation_2op_1t(
        M::AbstractHEOMLSMatrix,
        state::Union{QuantumObject,ADOs},
        τlist::AbstractVector,
        A::QuantumObject{Operator},
        B::QuantumObject{Operator};
        alg::AbstractODEAlgorithm = DP5(),
        reverse::Bool = false,
        matrix_form::Union{Val,Bool} = Val(false),
        progress_bar::Union{Val,Bool} = Val(true),
        kwargs...,
    )

Returns the two-time correlation function (with only one time coordinate ``\tau``) of two operators ``\hat{A}`` and ``\hat{B}`` : ``\left\langle \hat{A}(\tau) \hat{B}(0) \right\rangle`` for a given `state` using HEOM approach.

# Notes

- When `reverse=true`, the correlation function is calculated as ``\left\langle \hat{A}(t) \hat{B}(t + \tau) \right\rangle``.
"""
function QuantumToolbox.correlation_2op_1t(
        M::AbstractHEOMLSMatrix,
        state::T_state,
        τlist::AbstractVector,
        A::QuantumObject{Operator},
        B::QuantumObject{Operator};
        alg::AbstractODEAlgorithm = DP5(),
        reverse::Bool = false,
        progress_bar::Union{Val, Bool} = Val(true),
        kwargs...,
    ) where {T_state <: Union{QuantumObject, ADOs}}
    # there is only a single t = 0, so the underlying HEOMsolve_map ensemble has just one trajectory
    corr = correlation_2op_2t(M, state, [0], τlist, A, B; alg, reverse, ensemblealg = EnsembleSerial(), progress_bar, kwargs...)

    return corr[1, :] # 1 means tlist[1] = 0
end

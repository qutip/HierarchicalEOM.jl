@doc raw"""
    correlation_3op_2t(
        M::AbstractHEOMLSMatrix,
        state::Union{QuantumObject,ADOs},
        tlist::AbstractVector,
        τlist::AbstractVector,
        A::QuantumObject{Operator},
        B::QuantumObject{Operator},
        C::QuantumObject{Operator};
        kwargs...,
    )

Returns the two-times correlation function of three operators ``\hat{A}``, ``\hat{B}`` and ``\hat{C}``: ``\left\langle \hat{A}(t) \hat{B}(t + \tau) \hat{C}(t) \right\rangle`` for a given `state` using HEOM approach.
"""
function QuantumToolbox.correlation_3op_2t(
    M::AbstractHEOMLSMatrix,
    state::T_state,
    tlist::AbstractVector,
    τlist::AbstractVector,
    A::QuantumObject{Operator},
    B::QuantumObject{Operator},
    C::QuantumObject{Operator};
    kwargs...,
) where {T_state<:Union{QuantumObject,ADOs}}

    # check tlist and τlist
    QuantumToolbox._check_correlation_time_list(tlist)
    QuantumToolbox._check_correlation_time_list(τlist)

    AC = HEOMSuperOp(sprepost(C, A), M.parity, M)

    kwargs2 = merge((saveat = collect(tlist),), (; kwargs...))
    ados_t_list = HEOMsolve(M, state, tlist; kwargs2...).ados

    corr = map((t, ρt) -> HEOMsolve(M, AC * ρt, τlist .+ t, e_ops = [B]; kwargs...).expect[1, :], tlist, ados_t_list)
    return reduce(vcat, transpose.(corr))
end

@doc raw"""
    correlation_3op_1t(
        M::AbstractHEOMLSMatrix,
        state::Union{QuantumObject,ADOs},
        τlist::AbstractVector,
        A::QuantumObject{Operator},
        B::QuantumObject{Operator},
        C::QuantumObject{Operator};
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
    kwargs...,
) where {T_state<:Union{QuantumObject,ADOs}}
    corr = correlation_3op_2t(M, state, [0], τlist, A, B, C; kwargs...)
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
        reverse::Bool = false,
        kwargs...,
    )

Returns the two-times correlation function of two operators ``\hat{A}`` and ``\hat{B}`` : ``\left\langle \hat{A}(t + \tau) \hat{B}(t) \right\rangle`` for a given `state` using HEOM approach.

When `reverse=true`, the correlation function is calculated as ``\left\langle \hat{A}(t) \hat{B}(t + \tau) \right\rangle``.
"""
function QuantumToolbox.correlation_2op_2t(
    M::AbstractHEOMLSMatrix,
    state::T_state,
    tlist::AbstractVector,
    τlist::AbstractVector,
    A::QuantumObject{Operator},
    B::QuantumObject{Operator};
    reverse::Bool = false,
    kwargs...,
) where {T_state<:Union{QuantumObject,ADOs}}
    C = eye(prod(M.dimensions), dims = M.dimensions)

    if reverse
        corr = correlation_3op_2t(M, state, tlist, τlist, A, B, C; kwargs...)
    else
        corr = correlation_3op_2t(M, state, tlist, τlist, C, A, B; kwargs...)
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
        reverse::Bool = false,
        kwargs...,
    )

Returns the two-time correlation function (with only one time coordinate ``\tau``) of two operators ``\hat{A}`` and ``\hat{B}`` : ``\left\langle \hat{A}(\tau) \hat{B}(0) \right\rangle`` for a given `state` using HEOM approach.

When `reverse=true`, the correlation function is calculated as ``\left\langle \hat{A}(0) \hat{B}(\tau) \right\rangle``.
"""
function QuantumToolbox.correlation_2op_1t(
    M::AbstractHEOMLSMatrix,
    state::T_state,
    τlist::AbstractVector,
    A::QuantumObject{Operator},
    B::QuantumObject{Operator};
    reverse::Bool = false,
    kwargs...,
) where {T_state<:Union{QuantumObject,ADOs}}
    corr = correlation_2op_2t(M, state, [0], τlist, A, B; reverse = reverse, kwargs...)

    return corr[1, :] # 1 means tlist[1] = 0
end


function _check_correlation_time_list(tlist::AbstractVector)
    any(t -> t == 0, tlist) ||
        throw(ArgumentError("The time list for calculating correlation function must contain the element `0`"))
    all(>=(0), tlist) ||
        throw(ArgumentError("All the elements in the time list for calculating correlation function must be positive."))
    return nothing
end

@doc raw"""
    QuantumToolbox.correlation_3op_2t(
        M::AbstractHEOMLSMatrix,
        ρ0::T_state,
        tlist::AbstractVector,
        τlist::AbstractVector,
        A::QuantumObject{Operator},
        B::QuantumObject{Operator},
        C::QuantumObject{Operator};
        kwargs...,
    ) where {T_state<:Union{QuantumObject,ADOs}}

``\left\langle \hat{A}(t) \hat{B}(t + \tau) \hat{C}(t) \right\rangle``
"""
function QuantumToolbox.correlation_3op_2t(
    M::AbstractHEOMLSMatrix,
    ρ0::T_state,
    tlist::AbstractVector,
    τlist::AbstractVector,
    A::QuantumObject{Operator},
    B::QuantumObject{Operator},
    C::QuantumObject{Operator};
    kwargs...,
) where {T_state<:Union{QuantumObject,ADOs}}

    # check tlist and τlist
    _check_correlation_time_list(tlist)
    _check_correlation_time_list(τlist)

    ρ0 = (T_state <: QuantumObject) ? ADOs(ρ0, M.N, M.parity) : ρ0
    
    AC = HEOMSuperOp(sprepost(C,A), M.parity, M)
    
    kwargs2 = merge((saveat = collect(tlist),), (; kwargs...))
    ρt_list = HEOMsolve(M, ρ0, tlist; kwargs2...).ados

    corr = map((t, ρt) -> HEOMsolve(M, AC * ρt, τlist .+ t, e_ops = [B]; kwargs...).expect[1, :], tlist, ρt_list)
    return reduce(vcat, transpose.(corr))
end

@doc raw"""
    QuantumToolbox.correlation_3op_1t(
        M::AbstractHEOMLSMatrix,
        ρ0::T_state,
        τlist::AbstractVector,
        A::QuantumObject{Operator},
        B::QuantumObject{Operator},
        C::QuantumObject{Operator};
        kwargs...,
    ) where {T_state<:Union{QuantumObject,ADOs}}

``\left\langle \hat{A}(0) \hat{B}(\tau) \hat{C}(0) \right\rangle`` for a given initial state ``|\psi_0\rangle``
"""
function QuantumToolbox.correlation_3op_1t(
    M::AbstractHEOMLSMatrix,
    ρ0::T_state,
    τlist::AbstractVector,
    A::QuantumObject{Operator},
    B::QuantumObject{Operator},
    C::QuantumObject{Operator};
    kwargs...,
) where {T_state<:Union{QuantumObject,ADOs}}
    corr = correlation_3op_2t(M, ρ0, [0], τlist, A, B, C)
    return corr[1,:]
end

@doc raw"""
    QuantumToolbox.correlation_2op_2t(
        M::AbstractHEOMLSMatrix,
        ρ0::T_state,
        tlist::AbstractVector,
        τlist::AbstractVector,
        A::QuantumObject{Operator},
        B::QuantumObject{Operator};
        reverse::Bool = false,
        kwargs...,
    ) where {T_state<:Union{QuantumObject,ADOs}}

``\left\langle \hat{A}(t + \tau) \hat{B}(t) \right\rangle``
"""
function QuantumToolbox.correlation_2op_2t(
    M::AbstractHEOMLSMatrix,
    ρ0::T_state,
    tlist::AbstractVector,
    τlist::AbstractVector,
    A::QuantumObject{Operator},
    B::QuantumObject{Operator};
    reverse::Bool = false,
    kwargs...,
) where {T_state<:Union{QuantumObject,ADOs}}
    C = eye(prod(M.dimensions), dims = M.dimensions)
    
    if reverse
        corr = correlation_3op_2t(M, ρ0, tlist, τlist, A, B, C; kwargs...)
    else
        corr = correlation_3op_2t(M, ρ0, tlist, τlist, C, A, B; kwargs...)
    end

    return corr
end

@doc raw"""
    QuantumToolbox.correlation_2op_1t(
        M::AbstractHEOMLSMatrix,
        ρ0::T_state,
        τlist::AbstractVector,
        A::QuantumObject{Operator},
        B::QuantumObject{Operator};
        reverse::Bool = false,
        kwargs...,
    ) where {T_state<:Union{QuantumObject,ADOs}}

``\left\langle \hat{A}(\tau) \hat{B}(0) \right\rangle``

if `reverse` ``\left\langle \hat{A}(0) \hat{B}(\tau) \right\rangle``
"""
function QuantumToolbox.correlation_2op_1t(
    M::AbstractHEOMLSMatrix,
    ρ0::T_state,
    τlist::AbstractVector,
    A::QuantumObject{Operator},
    B::QuantumObject{Operator};
    reverse::Bool = false,
    kwargs...,
) where {T_state<:Union{QuantumObject,ADOs}}

    corr = correlation_2op_2t(M, ρ0, [0], τlist, A, B; reverse = reverse, kwargs...)

    return corr[1, :] # 1 means tlist[1] = 0
end

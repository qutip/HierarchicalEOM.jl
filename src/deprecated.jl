#=
This file gathers all the deprecated names (structures, functions, or variables) which will be removed in the future major release.

- Before the major release, the deprecated names will just throw errors when they are called.
- If the deprecated names were once exported, we will still export them here until next major release.
- If we decide to push a major release, cleanup this file.

Example:

export deprecated_foo

function deprecated_foo(args...; kwargs...)
    error("`deprecated_foo` has been deprecated and will be removed in next major release, please use `new_foo` instead.")
end
=#

export evolution, SteadyState, C

SteadyState(args...; kwargs...) = error("`SteadyState` has been deprecated, please use `steadystate` instead.")

evolution(args...; kwargs...) = error("`evolution` has been deprecated, please use `HEOMsolve` instead.")

C(args...; kwargs...) = error("`C` has been deprecated, please use `correlation_function` instead.")

_HEOMSuperOp_deprecated_method_error() = error(
    "Specifying `mul_basis = \"L\", \"R\", or \"LR\"` to `HEOMSuperOp` has been deprecated, try to construct system SuperOperator manually by using `spre`, `spost`, or `sprepost`.",
)

HEOMSuperOp(op, opParity::AbstractParity, dims, N::Int, mul_basis::AbstractString; Id_cache = I(N)) =
    _HEOMSuperOp_deprecated_method_error()
HEOMSuperOp(
    op,
    opParity::AbstractParity,
    refHEOMLS::AbstractHEOMLSMatrix,
    mul_basis::AbstractString;
    Id_cache = I(refHEOMLS.N),
) = HEOMSuperOp(op, opParity, refHEOMLS.dimensions, refHEOMLS.N, mul_basis; Id_cache = Id_cache)
HEOMSuperOp(op, opParity::AbstractParity, refADOs::ADOs, mul_basis::AbstractString; Id_cache = I(refADOs.N)) =
    HEOMSuperOp(op, opParity, refADOs.dimensions, refADOs.N, mul_basis; Id_cache = Id_cache)

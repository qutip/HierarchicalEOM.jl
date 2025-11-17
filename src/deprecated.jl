#=
This file gathers all the deprecated names (structures, functions, or variables) which will be removed in the future major release.

- Before the major release, the deprecated names will show warnings or just throw errors when they are called.
- If the deprecated names were once exported, we will still export them here until next major release.
- If we decide to push a major release, cleanup this file.

Example 1 [throw errors if the deprecation is fundamental, and there is no replacement for it]:
```
export deprecated_foo
deprecated_foo(args...; kwargs...) = error("`deprecated_foo` is deprecated and will be removed in next major release.")
```

Example 2 ["force" show warning and tell the users that there is a replacement for the deprecated function]:

```
export deprecated_foo
function deprecated_foo(args...; kwargs...)
    Base.depwarn("`deprecated_foo` is deprecated and will be removed in next major release, use `new_foo` instead.", :deprecated_foo, force = true)
    new_foo(args...; kwargs...)
end
```
=#

export SteadyState
function SteadyState(args...; kwargs...)
    Base.depwarn(
        "`SteadyState` is deprecated and will be removed in next major release, use `steadystate` instead.",
        :SteadyState,
        force = true,
    )
    return steadystate(args...; kwargs...)
end

export evolution
function evolution(args...; kwargs...)
    Base.depwarn(
        "`evolution` is deprecated and will be removed in next major release, use `HEOMsolve` instead.",
        :evolution,
        force = true,
    )
    return HEOMsolve(args...; kwargs...)
end

export C
function C(args...)
    Base.depwarn(
        "`C` is deprecated and will be removed in next major release, use `correlation_function` instead.",
        :C,
        force = true,
    )
    return correlation_function(args...)
end

_HEOMSuperOp_deprecated_method_error() = error(
    "Specifying `mul_basis = \"L\", \"R\", or \"LR\"` to `HEOMSuperOp` is deprecated, try to construct system SuperOperator manually by using `spre`, `spost`, or `sprepost`.",
)

HEOMSuperOp(op, opParity::AbstractParity, dims, N::Int, mul_basis::AbstractString) =
    _HEOMSuperOp_deprecated_method_error()
HEOMSuperOp(
    op,
    opParity::AbstractParity,
    refHEOMLS::AbstractHEOMLSMatrix,
    mul_basis::AbstractString
) = HEOMSuperOp(op, opParity, refHEOMLS.dimensions, refHEOMLS.N, mul_basis)
HEOMSuperOp(op, opParity::AbstractParity, refADOs::ADOs, mul_basis::AbstractString) =
    HEOMSuperOp(op, opParity, refADOs.dimensions, refADOs.N, mul_basis)

export AbstractParity
export OddParity, ODD
export EvenParity, EVEN

abstract type AbstractParity end

@doc raw"""
    struct OddParity <: AbstractParity
"""
struct OddParity <: AbstractParity end

@doc raw"""
    struct EvenParity <: AbstractParity
"""
struct EvenParity <: AbstractParity end

value(p::OddParity) = 1
value(p::EvenParity) = 0

Base.:(!)(p::OddParity) = EVEN
Base.:(!)(p::EvenParity) = ODD
Base.:(*)(p1::TP1, p2::TP2) where {TP1,TP2<:AbstractParity} = TP1 == TP2 ? EVEN : ODD

Base.show(io::IO, p::OddParity) = print(io, "odd-parity")

Base.show(io::IO, p::EvenParity) = print(io, "even-parity")

# Parity label
@doc raw"""
    const ODD  = OddParity()
Label of odd-parity
"""
const ODD = OddParity()

@doc raw"""
    const EVEN = EvenParity()
Label of even-parity
"""
const EVEN = EvenParity()

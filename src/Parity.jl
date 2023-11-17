abstract type AbstractParity end

@doc raw"""
    struct OddParity <: AbstractParity
"""
struct OddParity  <: AbstractParity end

@doc raw"""
    struct EvenParity <: AbstractParity
"""
struct EvenParity <: AbstractParity end

value(p::OddParity)  = 1
value(p::EvenParity) = 0

!(p::OddParity)  = EVEN
!(p::EvenParity) = ODD
*(p1::TP1, p2::TP2) where {TP1, TP2 <: AbstractParity} = TP1 == TP2 ? EVEN : ODD

function show(io::IO, p::OddParity)
    print(io, "odd-parity")
end

function show(io::IO, p::EvenParity)
    print(io, "even-parity")
end

# Parity label
@doc raw"""
    const ODD  = OddParity()
Label of odd-parity
"""
const ODD  = OddParity()

@doc raw"""
    const EVEN = EvenParity()
Label of even-parity
"""
const EVEN = EvenParity()
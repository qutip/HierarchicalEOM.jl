@testitem "Nvec equality and hashing" begin
    using SparseArrays
    using HierarchicalEOM: Nvec

    # `==` must compare the actual data, not the hashes. Defining equality via hash
    # equality breaks the `hash`/`==` contract that `Dict` relies on: a 64-bit hash
    # collision between two distinct ADOs would then be declared "equal", conflating
    # them into a single key and corrupting the hierarchy mapping (issue #298).

    # equal data => equal, and (necessarily) equal hash
    a = Nvec([0, 1, 2, 0, 3])
    b = Nvec([0, 1, 2, 0, 3])
    @test a == b
    @test hash(a) == hash(b)

    # different data (even a single neighbouring index, as in the coupling terms)
    # must compare unequal
    c = Nvec([0, 1, 2, 0, 3])
    d = Nvec([0, 1, 3, 0, 2])  # same level, different distribution
    @test c != d
    e = Nvec([0, 1, 2, 0, 4])  # differs by one excitation at a single index
    @test c != e

    # sanity: equality is reflexive and consistent with copy
    @test a == copy(a)

    # Nvecs must behave as distinct keys in a Dict (the core use in HierarchyDict)
    dict = Dict{Nvec,Int}()
    dict[c] = 1
    dict[d] = 2
    dict[e] = 3
    @test length(dict) == 3
    @test dict[Nvec([0, 1, 2, 0, 3])] == 1
end

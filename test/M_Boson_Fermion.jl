@testitem "M_Boson_Fermion" begin
    using SparseArrays
    using SciMLOperators

    # Test Boson-Fermion-type HEOM Liouvillian superoperator matrix
    λ = 0.1450
    W = 0.6464
    kT = 0.7414
    μ = 0.8787
    N = 3
    tierb = 2
    tierf = 2

    # System Hamiltonian
    Hsys = Qobj([
        0.6969 0.4364
        0.4364 0.3215
    ])

    # system-bath coupling operator
    Q = Qobj([
        0.1234 0.1357+0.2468im
        0.1357-0.2468im 0.5678
    ])
    Bbath = Boson_DrudeLorentz_Pade(Q, λ, W, kT, N)
    Fbath = Fermion_Lorentz_Pade(Q, λ, μ, W, kT, N)

    # jump operator
    J = Qobj([0 0.1450-0.7414im; 0.1450+0.7414im 0])

    L = M_Boson_Fermion(Hsys, tierb, tierf, Bbath, Fbath; verbose = true) # also test verbosity
    L_combine = M_Boson_Fermion(Hsys, tierb, tierf, Bbath, Fbath; verbose = false, assemble = Val(:combine)) # test combined tensor assembly with assemble = Val(:combine)
    L_lazy = M_Boson_Fermion(Hsys, tierb, tierf, Bbath, Fbath; verbose = false, assemble = Val(:none)) # test lazy tensor assembly with assemble = Val(:none)
    @test show(devnull, MIME("text/plain"), L) === nothing
    @test size(L) == (2220, 2220)
    @test L.N == 555
    @test nnz(L.data.A) == nnz(L(0).data.A) == nnz(concretize(L_lazy.data)) == 43368
    @test L.data isa SciMLOperators.MatrixOperator
    @test L_combine.data isa SciMLOperators.AddedOperator
    @test L_lazy.data isa SciMLOperators.AddedOperator
    L = addBosonDissipator(L, J)
    @test nnz(L.data.A) == nnz(L(0).data.A) == 45590
    @test isconstant(L)
    @test iscached(L)
    ados = steadystate(L; verbose = false)
    @test ados.dims == L.dims
    @test length(ados) == L.N
    @test eltype(L) == eltype(ados)
    ρ0 = ados[1]
    @test getRho(ados) == ρ0
    ρ1 = Qobj([
        0.496709+4.88415e-12im -0.00324048+0.00286376im
        -0.00324048-0.00286376im 0.503291-4.91136e-12im
    ])
    @test ρ0 ≈ ρ1 atol=1e-6

    L = M_Boson_Fermion(Hsys, tierb, tierf, [Bbath, Bbath], Fbath; verbose = false)
    @test size(L) == (6660, 6660)
    @test L.N == 1665
    @test nnz(L.data.A) == nnz(L(0).data.A) == 139210
    L = addFermionDissipator(L, J)
    @test nnz(L.data.A) == nnz(L(0).data.A) == 145872
    ados = steadystate(L; verbose = false)
    @test ados.dims == L.dims
    @test length(ados) == L.N
    ρ0 = ados[1]
    @test getRho(ados) == ρ0
    ρ1 = Qobj([
        0.493774+6.27624e-13im -0.00536526+0.00651746im
        -0.00536526-0.00651746im 0.506226-6.15855e-13im
    ])
    @test ρ0 ≈ ρ1 atol=1e-6

    L = M_Boson_Fermion(Hsys, tierb, tierf, Bbath, [Fbath, Fbath]; verbose = false)
    @test size(L) == (8220, 8220)
    @test L.N == 2055
    @test nnz(L.data.A) == nnz(L(0).data.A) == 167108
    L = addBosonDissipator(L, J)
    @test nnz(L.data.A) == nnz(L(0).data.A) == 175330
    ados = steadystate(L; verbose = false)
    @test ados.dims == L.dims
    @test length(ados) == L.N
    ρ0 = ados[1]
    @test getRho(ados) == ρ0
    ρ1 = Qobj([
        0.496468-4.32253e-12im -0.00341484+0.00316445im
        -0.00341484-0.00316445im 0.503532+4.32574e-12im
    ])
    @test ρ0 ≈ ρ1 atol=1e-6

    ## check exceptions
    @test_throws ErrorException ados[L.N+1]
    @test_throws ErrorException M_Boson_Fermion(Qobj([0, 0]), tierb, tierf, Bbath, Fbath; verbose = false)
end

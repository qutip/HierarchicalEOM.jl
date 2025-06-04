@testitem "M_Boson" begin
    using SparseArrays

    # Test Boson-type HEOM Liouvillian superoperator matrix
    λ = 0.1450
    W = 0.6464
    kT = 0.7414
    μ = 0.8787
    N = 5
    tier = 3

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

    # jump operator
    J = Qobj([0 0.1450-0.7414im; 0.1450+0.7414im 0])

    L = M_Boson(Hsys, tier, Bbath; verbose = false)
    @test show(devnull, MIME("text/plain"), L) === nothing
    @test size(L) == (336, 336)
    @test L.N == 84
    @test nnz(L.data.A) == nnz(L(0)) == 4422
    L = addBosonDissipator(L, J)
    @test nnz(L.data.A) == nnz(L(0)) == 4760
    ados = steadystate(L; verbose = false)
    @test ados.dims == L.dims
    @test length(ados) == L.N
    @test eltype(L) == eltype(ados)
    ρ0 = ados[1]
    @test getRho(ados) == ρ0
    ρ1 = Qobj(
        [
            0.4969521584882579-2.27831302340618e-13im -0.0030829715611090133+0.002534368458048467im
            -0.0030829715591718203-0.0025343684616701547im 0.5030478415140676+2.3661885315257474e-13im
        ],
    )
    @test ρ0 ≈ ρ1

    L = M_Boson(Hsys, tier, [Bbath, Bbath]; verbose = false)
    @test size(L) == (1820, 1820)
    @test L.N == 455
    @test nnz(L.data.A) == nnz(L(0)) == 27662
    L = addBosonDissipator(L, J)
    @test nnz(L.data.A) == nnz(L(0)) == 29484
    ados = steadystate(L; verbose = false)
    @test ados.dims == L.dims
    @test length(ados) == L.N
    ρ0 = ados[1]
    @test getRho(ados) == ρ0
    ρ1 = Qobj(
        [
            0.49406682844513267+9.89558173111355e-13im -0.005261234545120281+0.0059968903987593im
            -0.005261234550122085-0.005996890386139547im 0.5059331715578721-9.413847493320824e-13im
        ],
    )
    @test ρ0 ≈ ρ1

    ## check exceptions
    @test_throws BoundsError L[1, 1821]
    @test_throws BoundsError L[1821, 336]
    @test_throws ErrorException ados[L.N+1]
    @test_throws ErrorException M_Boson(Qobj([0, 0]), tier, Bbath; verbose = false)
end

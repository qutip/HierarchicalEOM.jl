@testitem "M_Fermion" begin
    using SparseArrays
    using SciMLOperators

    # Test Fermion-type HEOM Liouvillian superoperator matrix
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
    Fbath = Fermion_Lorentz_Pade(Q, λ, μ, W, kT, N)

    # jump operator
    J = Qobj([0 0.1450-0.7414im; 0.1450+0.7414im 0])

    L = M_Fermion(Hsys, tier, Fbath; verbose = true) # also test verbosity
    L_lazy = M_Fermion(Hsys, tier, Fbath; verbose = false, concretize = Val(false))
    @test show(devnull, MIME("text/plain"), L) === nothing
    @test size(L) == (1196, 1196)
    @test L.N == 299
    @test nnz(L.data.A) == nnz(L(0)) == nnz(concretize(L_lazy)(0)) == 21318
    @test L_lazy.data isa SciMLOperators.AddOperator
    @test length(L_lazy.data.ops) == 8 * 1 + 2 # 8 ops per fermion bath + 2 free terms
    L = addFermionDissipator(L, J)
    @test nnz(L.data.A) == nnz(L(0)) == 22516
    @test isconstant(L)
    @test iscached(L)
    ados = steadystate(L; verbose = false)
    @test ados.dims == L.dims
    @test length(ados) == L.N
    @test eltype(L) == eltype(ados)
    ρ0 = ados[1]
    @test getRho(ados) == ρ0
    ρ1 = Qobj(
        [
            0.49971864340781574+1.5063528845574697e-11im -0.00025004129095353573+0.00028356932981729176im
            -0.0002500413218393161-0.0002835693203755187im 0.5002813565929579-1.506436545359778e-11im
        ],
    )
    @test ρ0 ≈ ρ1 atol=1e-6

    L = M_Fermion(Hsys, tier, Fbath; threshold = 1e-8, verbose = false)
    L = addFermionDissipator(L, J)
    @test size(L) == (148, 148)
    @test L.N == 37
    @test nnz(L.data.A) == nnz(L(0)) == 2054
    ados = steadystate(L; verbose = false)
    ρ2 = ados[1]
    @test ρ0 ≈ ρ2 atol=1e-6

    L = M_Fermion(Hsys, tier, [Fbath, Fbath]; verbose = false)
    @test size(L) == (9300, 9300)
    @test L.N == 2325
    @test nnz(L.data.A) == nnz(L(0)) == 174338
    L = addFermionDissipator(L, J)
    @test nnz(L.data.A) == nnz(L(0)) == 183640
    ados = steadystate(L; verbose = false)
    @test ados.dims == L.dims
    @test length(ados) == L.N
    ρ0 = ados[1]
    @test getRho(ados) == ρ0
    ρ1 = Qobj(
        [
            0.4994229368103249+2.6656157051929377e-12im -0.0005219753638749304+0.0005685093274121244im
            -0.0005219753958601764-0.0005685092413099392im 0.5005770631903512-2.6376966158390854e-12im
        ],
    )
    @test ρ0 ≈ ρ1

    ## check exceptions
    @test_throws BoundsError L[1, 9301]
    @test_throws BoundsError L[9301, 9300]
    @test_throws ErrorException ados[L.N+1]
    @test_throws ErrorException M_Fermion(Qobj([0, 0]), tier, Fbath; verbose = false)
end

@testitem "M_Fermion" begin
    using SparseArrays
    using SciMLOperators

    # Test Fermion-type HEOM Liouvillian superoperator matrix
    λ = 0.145
    W = 0.6464
    kT = 0.7414
    μ = 0.8787
    N = 5
    tier = 3

    # System Hamiltonian
    Hsys = Qobj(
        sparse(
            [
                0.6969 0.4364
                0.4364 0.3215
            ]
        )
    )

    # system-bath coupling operator
    Fbath = Fermion_Lorentz_Pade(destroy(2), λ, μ, W, kT, N)

    # jump operator
    J = Qobj(
        sparse(
            [
                0 0.145 - 0.7414im
                0.145 + 0.7414im 0
            ]
        )
    )

    L = M_Fermion(Hsys, tier, Fbath; verbose = true) # also test verbosity
    L_combine = M_Fermion(Hsys, tier, Fbath; verbose = false, assemble = Val(:combine))
    L_lazy = M_Fermion(Hsys, tier, Fbath; verbose = false, assemble = Val(:none))
    L_combine_cached = cache_operator(L_combine, similar(zeros(eltype(L_combine), size(L_combine, 1))))
    @test show(devnull, MIME("text/plain"), L) === nothing
    @test size(L) == (1196, 1196)
    @test L.N == 299
    @test nnz(L.data.A) == nnz(L(0).data.A) == nnz(concretize(L_lazy)(0).data.A) == 10018
    @test L(0).data.A == concretize(L_combine).data.A
    @test L.data isa SciMLOperators.MatrixOperator
    @test issparse(L.data.A) # check if it's a sparse matrix
    @test L_combine.data isa SciMLOperators.AddedOperator
    @test L_lazy.data isa SciMLOperators.AddedOperator
    @test length(L_combine.data.ops) == 4 * 1 + 2 # 4 ops per fermion bath + 2 free terms
    @test length(L_lazy.data.ops) == 8 * 1 + 2 # 8 ops per fermion bath + 2 free terms
    L = addFermionDissipator(L, J)
    @test nnz(L.data.A) == nnz(L(0).data.A) == 11216
    @test isconstant(L)
    @test iscached(L)
    @test iscached(L_combine_cached)
    ados = steadystate(L; verbose = false)
    @test ados.dims.to == L.dims.to
    @test length(ados) == L.N
    @test eltype(L) == eltype(ados)
    ρ0 = ados[1]
    @test getRho(ados) == ρ0
    ρ1 = Qobj(
        [
            0.5014373576999699 - 1.3892567817792988e-17im  -0.0015662233262583364 - 0.009015267133592107im
            -0.0015662233262584758 + 0.009015267133592103im      0.4985626423000295 - 6.95905692573837e-18im
        ],
    )
    @test ρ0 ≈ ρ1 atol = 1.0e-6

    L = M_Fermion(Hsys, tier, Fbath; threshold = 1.0e-8, verbose = false)
    L = addFermionDissipator(L, J)
    @test size(L) == (148, 148)
    @test L.N == 37
    @test nnz(L.data.A) == nnz(L(0).data.A) == 1112
    ados = steadystate(L; verbose = false)
    ρ2 = ados[1]
    @test ρ0 ≈ ρ2 atol = 1.0e-6

    L = M_Fermion(Hsys, tier, [Fbath, Fbath]; verbose = false)
    @test size(L) == (9300, 9300)
    @test L.N == 2325
    @test nnz(L.data.A) == nnz(L(0).data.A) == 81082
    L = addFermionDissipator(L, J)
    @test nnz(L.data.A) == nnz(L(0).data.A) == 90384
    ados = steadystate(L; verbose = false)
    @test ados.dims.to == L.dims.to
    @test length(ados) == L.N
    ρ0 = ados[1]
    @test getRho(ados) == ρ0
    ρ1 = Qobj(
        [
            0.5024605153620629 - 3.452236695213822e-17im  -0.003127033937084677 - 0.017184962302082486im
            -0.0031270339370842257 + 0.01718496230208245im        0.4975394846379351 + 3.9008249021550266e-17im
        ],
    )
    @test ρ0 ≈ ρ1

    ## check exceptions
    @test_throws BoundsError L[1, 9301]
    @test_throws BoundsError L[9301, 9300]
    @test_throws ErrorException ados[L.N + 1]
    @test_throws ErrorException M_Fermion(Qobj([0, 0]), tier, Fbath; verbose = false)
end

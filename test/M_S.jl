@testitem "M_S" begin
    using SparseArrays

    # Test Schrodinger type HEOM Liouvillian superoperator matrix
    t = 10
    Hsys = sigmax()
    L = M_S(Hsys; verbose = true) # also test verbosity
    L_super = M_S(liouvillian(Hsys); verbose = false) # test if input is already Liouvillian
    ψ0 = basis(2, 0)
    @test show(devnull, MIME("text/plain"), L) === nothing
    @test size(L) == size(L_super) == (4, 4)
    @test L.N == L_super.N == 1
    @test nnz(L.data.A) == nnz(L_super.data.A) == nnz(L(0).data.A) == nnz(L_super(0).data.A) == 8
    @test L.data == L_super.data
    @test isconstant(L)
    @test iscached(L)
    ados_list = heomsolve(L, ψ0, 0:1:t; reltol = 1.0e-8, abstol = 1.0e-10, progress_bar = Val(false)).ados
    ados = ados_list[end]
    @test ados.dims == L.dims
    @test length(ados) == L.N
    @test eltype(L) == eltype(ados)
    ρ1 = ados[1]
    @test getRho(ados) == ρ1
    ρ2 = Qobj(
        [
            cos(t)^2 0.5im * sin(2 * t)
            -0.5im * sin(2 * t) sin(t)^2
        ]
    )
    @test ρ1 ≈ ρ2 atol = 1.0e-6

    L = addBosonDissipator(L, √(0.01) * sigmaz())
    L = addFermionDissipator(L, √(0.01) * sigmaz())
    @test nnz(L.data.A) == nnz(L(0).data.A) == 10
    ados_list = heomsolve(L, ψ0, 0:0.5:t; reltol = 1.0e-8, abstol = 1.0e-10, progress_bar = Val(false)).ados
    ados = ados_list[end]
    ρ3 = ados[1]
    @test getRho(ados) == ρ3
    ρ4 = Qobj(
        [
            0.6711641119639493 + 0.0im 0.3735796014268062im
            -0.3735796014268062im 0.32883588803605024
        ]
    )
    @test ρ3 ≈ ρ4

    ## check exceptions
    @test_throws BoundsError L[1, 5]
    @test_throws BoundsError L[5, 2]
    @test_throws ErrorException ados[L.N + 1]
    @test_throws ErrorException M_S(Qobj([0, 0]); verbose = false)
end

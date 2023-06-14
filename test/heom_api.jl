# Test Schrodinger type HEOM Liouvillian superoperator matrix
@testset "M_S" begin
    t = 10
    Hsys = [0 1; 1 0]
    L = M_S(Hsys; verbose=false)
    @test show(devnull, MIME("text/plain"), L) == nothing
    @test size(L) == (4, 4)
    @test L.N  == 1
    @test nnz(L.data) == 8
    ados_list = evolution(L, [1 0; 0 0], 0:1:t; reltol=1e-8, abstol=1e-10, verbose=false)
    ados = ados_list[end]
    @test ados.dim == L.dim
    @test length(ados) == L.N
    ρ0 = ados[1]
    @test getRho(ados) == ρ0
    ρ1 = [
              cos(t)^2 -0.5im*sin(2*t);
        0.5im*sin(2*t)        sin(t)^2
    ]
    @test _is_Matrix_approx(ρ0, ρ1)

    L =   addBosonDissipator(L, √(0.01) * [1 0; 0 -1])
    L = addFermionDissipator(L, √(0.01) * [1 0; 0 -1])
    @test nnz(L.data) == 10
    ados_list = evolution(L, [1 0; 0 0], 0:0.5:t; reltol=1e-8, abstol=1e-10, verbose=false)
    ados = ados_list[end]
    ρ0 = ados[1]
    @test getRho(ados) == ρ0
    ρ1 = [
        0.6711641119639493+0.0im 0.3735796014268062im;
        -0.3735796014268062im    0.32883588803605024
    ]
    @test _is_Matrix_approx(ρ0, ρ1)

    ## check exceptions
    @test_throws BoundsError L[1, 5]
    @test_throws BoundsError L[1:5, 2]
    @test_throws ErrorException ados[L.N + 1]
    @test_throws ErrorException M_S(Hsys, :wrong; verbose=false)
    @test_throws ErrorException @test_warn "HEOM doesn't support matrix type : Vector{Int64}" M_S([0, 0]; verbose=false)
end

λ = 0.1450
W = 0.6464
T = 0.7414
μ = 0.8787
N = 5
tier = 3

# System Hamiltonian
Hsys = [
    0.6969 0.4364;
    0.4364 0.3215
]

# system-bath coupling operator
Q = [
               0.1234 0.1357 + 0.2468im; 
    0.1357 - 0.2468im            0.5678
]
Bbath = Boson_DrudeLorentz_Pade(Q, λ, W, T, N)
Fbath = Fermion_Lorentz_Pade(Q, λ, μ, W, T, N)

# jump operator
J = [0 0.1450 - 0.7414im; 0.1450 + 0.7414im 0]

# Test Boson-type HEOM Liouvillian superoperator matrix
@testset "M_Boson" begin
    L = M_Boson(Hsys, tier, Bbath; verbose=false)
    @test show(devnull, MIME("text/plain"), L) == nothing
    @test size(L) == (336, 336)
    @test L.N  == 84
    @test nnz(L.data) == 4422
    L = addBosonDissipator(L, J)
    @test nnz(L.data) == 4760
    ados = SteadyState(L, [0.64 0; 0 0.36]; verbose=false)
    @test ados.dim == L.dim
    @test length(ados) == L.N
    ρ0 = ados[1]
    @test getRho(ados) == ρ0
    ρ1 = [
        0.4969521584882579 - 2.27831302340618e-13im -0.0030829715611090133 + 0.002534368458048467im; 
        -0.0030829715591718203 - 0.0025343684616701547im 0.5030478415140676 + 2.3661885315257474e-13im
    ]
    @test _is_Matrix_approx(ρ0, ρ1)

    L = M_Boson(Hsys, tier, [Bbath, Bbath]; verbose=false)
    @test size(L) == (1820, 1820)
    @test L.N  == 455
    @test nnz(L.data) == 27662
    L = addBosonDissipator(L, J)
    @test nnz(L.data) == 29484
    ados = SteadyState(L, [0.64 0; 0 0.36]; verbose=false)
    @test ados.dim == L.dim
    @test length(ados) == L.N
    ρ0 = ados[1]
    @test getRho(ados) == ρ0
    ρ1 = [
        0.49406682844513267 + 9.89558173111355e-13im  -0.005261234545120281 + 0.0059968903987593im;
        -0.005261234550122085 - 0.005996890386139547im      0.5059331715578721 - 9.413847493320824e-13im
    ]
    @test _is_Matrix_approx(ρ0, ρ1)

    ## check exceptions
    @test_throws BoundsError L[1, 1821]
    @test_throws BoundsError L[1:1821, 336]
    @test_throws ErrorException ados[L.N + 1]
    @test_throws ErrorException M_Boson(Hsys, tier, Bbath, :wrong; verbose=false)
    @test_throws ErrorException @test_warn "HEOM doesn't support matrix type : Vector{Int64}" M_Boson([0, 0], tier, Bbath; verbose=false)
end

# Test Fermion-type HEOM Liouvillian superoperator matrix
@testset "M_Fermion" begin
    L = M_Fermion(Hsys, tier, Fbath; verbose=false)
    @test show(devnull, MIME("text/plain"), L) == nothing
    @test size(L) == (1196, 1196)
    @test L.N  == 299
    @test nnz(L.data) == 21318
    L = addFermionDissipator(L, J)
    @test nnz(L.data) == 22516
    ados = SteadyState(L, [0.64 0; 0 0.36]; verbose=false)
    @test ados.dim == L.dim
    @test length(ados) == L.N
    ρ0 = ados[1]
    @test getRho(ados) == ρ0
    ρ1 = [
        0.49971864340781574 + 1.5063528845574697e-11im  -0.00025004129095353573 + 0.00028356932981729176im;
        -0.0002500413218393161 - 0.0002835693203755187im      0.5002813565929579 - 1.506436545359778e-11im
    ]
    @test _is_Matrix_approx(ρ0, ρ1)

    L = M_Fermion(Hsys, tier, Fbath; threshold = 1e-8, verbose=false)
    L = addFermionDissipator(L, J)
    @test size(L) == (148, 148)
    @test L.N  == 37
    @test nnz(L.data) == 2054
    ados = SteadyState(L, [0.64 0; 0 0.36]; verbose=false)
    ρ2 = ados[1]
    @test _is_Matrix_approx(ρ0, ρ2)

    L = M_Fermion(Hsys, tier, [Fbath, Fbath]; verbose=false)
    @test size(L) == (9300, 9300)
    @test L.N  == 2325
    @test nnz(L.data) == 174338
    L = addFermionDissipator(L, J)
    @test nnz(L.data) == 183640
    ados = SteadyState(L, [0.64 0; 0 0.36]; verbose=false)
    @test ados.dim == L.dim
    @test length(ados) == L.N
    ρ0 = ados[1]
    @test getRho(ados) == ρ0
    ρ1 = [
        0.4994229368103249 + 2.6656157051929377e-12im  -0.0005219753638749304 + 0.0005685093274121244im;
        -0.0005219753958601764 - 0.0005685092413099392im      0.5005770631903512 - 2.6376966158390854e-12im
    ]
    @test _is_Matrix_approx(ρ0, ρ1)

    ## check exceptions
    @test_throws BoundsError L[1, 9301]
    @test_throws BoundsError L[1:9301, 9300]
    @test_throws ErrorException ados[L.N + 1]
    @test_throws ErrorException M_Fermion(Hsys, tier, Fbath, :wrong; verbose=false)
    @test_throws ErrorException @test_warn "HEOM doesn't support matrix type : Vector{Int64}" M_Fermion([0, 0], tier, Fbath; verbose=false)
end

# Test Boson-Fermion-type HEOM Liouvillian superoperator matrix
@testset "M_Boson_Fermion" begin
    # re-define the bath (make the matrix smaller)
    λ = 0.1450
    W = 0.6464
    T = 0.7414
    μ = 0.8787
    N = 3
    tierb = 2
    tierf = 2

    Bbath = Boson_DrudeLorentz_Pade(Q, λ, W, T, N)
    Fbath = Fermion_Lorentz_Pade(Q, λ, μ, W, T, N)

    L = M_Boson_Fermion(Hsys, tierb, tierf, Bbath, Fbath; verbose=false)
    @test show(devnull, MIME("text/plain"), L) == nothing
    @test size(L) == (2220, 2220)
    @test L.N  == 555
    @test nnz(L.data) == 43368
    L = addBosonDissipator(L, J)
    @test nnz(L.data) == 45590
    ados = SteadyState(L, [0.64 0; 0 0.36]; verbose=false)
    @test ados.dim == L.dim
    @test length(ados) == L.N
    ρ0 = ados[1]
    @test getRho(ados) == ρ0
    ρ1 = [
        0.496709+4.88415e-12im  -0.00324048+0.00286376im;
    -0.00324048-0.00286376im      0.503291-4.91136e-12im
    ]
    @test _is_Matrix_approx(ρ0, ρ1)

    L = M_Boson_Fermion(Hsys, tierb, tierf, [Bbath, Bbath], Fbath; verbose=false)
    @test size(L) == (6660, 6660)
    @test L.N  == 1665
    @test nnz(L.data) == 139210
    L = addFermionDissipator(L, J)
    @test nnz(L.data) == 145872
    ados = SteadyState(L, [0.64 0; 0 0.36]; verbose=false)
    @test ados.dim == L.dim
    @test length(ados) == L.N
    ρ0 = ados[1]
    @test getRho(ados) == ρ0
    ρ1 = [
        0.493774+6.27624e-13im  -0.00536526+0.00651746im;
        -0.00536526-0.00651746im      0.506226-6.15855e-13im
    ]
    @test _is_Matrix_approx(ρ0, ρ1)

    L = M_Boson_Fermion(Hsys, tierb, tierf, Bbath, [Fbath, Fbath]; verbose=false)
    @test size(L) == (8220, 8220)
    @test L.N  == 2055
    @test nnz(L.data) == 167108
    L = addBosonDissipator(L, J)
    @test nnz(L.data) == 175330
    ados = SteadyState(L, [0.64 0; 0 0.36]; verbose=false)
    @test ados.dim == L.dim
    @test length(ados) == L.N
    ρ0 = ados[1]
    @test getRho(ados) == ρ0
    ρ1 = [
        0.496468-4.32253e-12im  -0.00341484+0.00316445im;
    -0.00341484-0.00316445im      0.503532+4.32574e-12im
    ]
    @test _is_Matrix_approx(ρ0, ρ1)

    ## check exceptions
    @test_throws ErrorException ados[L.N + 1]
    @test_throws ErrorException M_Boson_Fermion(Hsys, tierb, tierf, Bbath, Fbath, :wrong; verbose=false)
    @test_throws ErrorException @test_warn "HEOM doesn't support matrix type : Vector{Int64}" M_Boson_Fermion([0, 0], tierb, tierf, Bbath, Fbath; verbose=false)
end

@testset "Hierarchy Dictionary" begin
    Btier = 2
    Ftier = 2
    Nb = 3
    Nf = 3
    threshold = 1e-5

    Γ = 0.0025
    Dα  = 30
    Λ   = 0.0025
    ωcα = 0.2
    μL =  0.5
    μR = -0.5
    T = 0.025
    
    Hsys = [
        0   0     0     0;
        0 0.2     0     0;
        0   0 0.208  0.04;
        0   0  0.04 0.408
    ]

    ρ0 = [
        1 0 0 0;
        0 0 0 0;
        0 0 0 0;
        0 0 0 0
    ]

    cop = [
        0 1 0 0;
        1 0 0 0;
        0 0 0 1;
        0 0 1 0
    ]

    dop = [
        0 0 1 0;
        0 0 0 1;
        0 0 0 0;
        0 0 0 0
    ]

    bbath = Boson_DrudeLorentz_Matsubara(cop, Λ, ωcα, T, Nb)
    fbath = [
        Fermion_Lorentz_Pade(dop, Γ, μL, Dα, T, Nf),
        Fermion_Lorentz_Pade(dop, Γ, μR, Dα, T, Nf)
    ]

    L = M_Boson_Fermion(Hsys, Btier, Ftier, bbath, fbath; threshold = threshold, verbose=false)
    L = addTerminator(L, bbath)

    @test size(L) == (1696, 1696)
    @test L.N  == 106
    @test nnz(L.data) == 13536

    ados = SteadyState(L; verbose=false)
    @test Ic(ados, L, 1) ≈ 0.2883390125832726
    nvec_b, nvec_f = L.hierarchy.idx2nvec[1]
    @test_throws ErrorException getIndexEnsemble(nvec_f, L.hierarchy.bosonPtr)
    @test_throws ErrorException getIndexEnsemble(nvec_b, L.hierarchy.fermionPtr)
end

@testset "Auxiliary density operators" begin
    ados_b  = ADOs(spzeros(Int64, 20), 5)
    ados_f  = ADOs(spzeros(Int64,  8), 2)
    ados_bf = ADOs(spzeros(Int64, 40), 10)
    @test show(devnull, MIME("text/plain"), ados_b)  == nothing
    @test show(devnull, MIME("text/plain"), ados_f)  == nothing
    @test show(devnull, MIME("text/plain"), ados_bf) == nothing

    ρ_b = ados_b[:]
    # check iteration
    for (i, ado) in enumerate(ados_b)
       @test ρ_b[i] == ado
    end

    # expections for expect
    ados_wrong  = ADOs(spzeros(Int64, 18), 2)
    @test_throws ErrorException("The dimension of `op` is not consistent with `ados`.") @test_warn "The size of input matrix should be: (2, 2)." Expect([0 0 0; 0 0 0; 0 0 0], ados_f)
    @test_throws ErrorException("The dimension of the elements in `ados_list` should be consistent.") Expect([0 0; 0 0], [ados_b, ados_wrong])
    @test_throws ErrorException("The dimension of `op` is not consistent with the elements in `ados_list`.") @test_warn "The size of input matrix should be: (2, 2)." Expect([0 0 0; 0 0 0; 0 0 0], [ados_b, ados_f])
end
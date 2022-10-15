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

# Test Boson-type Heom liouvillian superoperator matrix
@testset "M_Boson" begin
    L = M_Boson(Hsys, tier, Bbath; verbose=false)
    @test size(L) == (336, 336)
    @test L.N  == 84
    @test L.Nb == 84
    @test L.Nf == 0
    @test nnz(L.data) == 4422
    addDissipator!(L, J)
    @test nnz(L.data) == 4760
    ados = SteadyState(L, [0.64 0; 0 0.36]; verbose=false)
    @test ados.dim == L.dim
    @test ados.Nb  == L.Nb
    @test ados.Nf  == L.Nf
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
    @test L.Nb == 455
    @test L.Nf == 0
    @test nnz(L.data) == 27662
    addDissipator!(L, J)
    @test nnz(L.data) == 29484
    ados = SteadyState(L, [0.64 0; 0 0.36]; verbose=false)
    @test ados.dim == L.dim
    @test ados.Nb  == L.Nb
    @test ados.Nf  == L.Nf
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
    @test_throws ErrorException ados[1, 1]
    @test_throws ErrorException getADO(ados, 1, 1)
    @test_throws ErrorException ados[L.N + 1]
    @test_throws ErrorException @test_warn "Heom doesn't support matrix type : Vector{Int64}" M_Boson([0, 0], tier, Bbath; verbose=false)
end

# Test Fermion-type Heom liouvillian superoperator matrix
@testset "M_Fermion" begin
    L = M_Fermion(Hsys, tier, Fbath; verbose=false)
    @test size(L) == (1196, 1196)
    @test L.N  == 299
    @test L.Nb == 0
    @test L.Nf == 299
    @test nnz(L.data) == 21318
    addDissipator!(L, J)
    @test nnz(L.data) == 22516
    ados = SteadyState(L, [0.64 0; 0 0.36]; verbose=false)
    @test ados.dim == L.dim
    @test ados.Nb  == L.Nb
    @test ados.Nf  == L.Nf
    @test length(ados) == L.N
    ρ0 = ados[1]
    @test getRho(ados) == ρ0
    ρ1 = [
        0.49971864340781574 + 1.5063528845574697e-11im  -0.00025004129095353573 + 0.00028356932981729176im;
        -0.0002500413218393161 - 0.0002835693203755187im      0.5002813565929579 - 1.506436545359778e-11im
    ]
    @test _is_Matrix_approx(ρ0, ρ1)

    L = M_Fermion(Hsys, tier, [Fbath, Fbath]; verbose=false)
    @test size(L) == (9300, 9300)
    @test L.N  == 2325
    @test L.Nb == 0
    @test L.Nf == 2325
    @test nnz(L.data) == 174338
    addDissipator!(L, J)
    @test nnz(L.data) == 183640
    ados = SteadyState(L, [0.64 0; 0 0.36]; verbose=false)
    @test ados.dim == L.dim
    @test ados.Nb  == L.Nb
    @test ados.Nf  == L.Nf
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
    @test_throws ErrorException ados[1, 1]
    @test_throws ErrorException getADO(ados, 1, 1)
    @test_throws ErrorException ados[L.N + 1]
    @test_throws ErrorException @test_warn "Heom doesn't support matrix type : Vector{Int64}" M_Fermion([0, 0], tier, Fbath; verbose=false)
end

# Test Boson-Fermion-type Heom liouvillian superoperator matrix
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
    @test size(L) == (2220, 2220)
    @test L.N  == 555
    @test L.Nb == 15
    @test L.Nf == 37
    @test nnz(L.data) == 43368
    addDissipator!(L, J)
    @test nnz(L.data) == 45590
    ados = SteadyState(L, [0.64 0; 0 0.36]; verbose=false)
    @test ados.dim == L.dim
    @test ados.Nb  == L.Nb
    @test ados.Nf  == L.Nf
    @test length(ados) == L.N
    ρ0 = ados[1, 1]
    @test getRho(ados) == ρ0
    ρ1 = [
        0.49693353824300623 - 1.1724586594620817e-7im  -0.0030854297725558563 + 0.0025734495019103824im;
        -0.003085174459117995 - 0.0025736446811396954im     0.5030664617592213 + 1.1724583626773739e-7im
    ]
    @test _is_Matrix_approx(ρ0, ρ1)

    L = M_Boson_Fermion(Hsys, tierb, tierf, [Bbath, Bbath], Fbath; verbose=false)
    @test size(L) == (6660, 6660)
    @test L.N  == 1665
    @test L.Nb == 45
    @test L.Nf == 37
    @test nnz(L.data) == 139210
    addDissipator!(L, J)
    @test nnz(L.data) == 145872
    ados = SteadyState(L, [0.64 0; 0 0.36]; verbose=false)
    @test ados.dim == L.dim
    @test ados.Nb  == L.Nb
    @test ados.Nf  == L.Nf
    @test length(ados) == L.N
    ρ0 = ados[1, 1]
    @test getRho(ados) == ρ0
    ρ1 = [
        0.49394119485917903 - 1.7614315992266511e-7im  -0.005274541015933129 + 0.006264966491795586im;
        -0.0052740922161028814 - 0.006265214322466733im     0.5060588051249636 + 1.761430060529106e-7im
    ]
    @test _is_Matrix_approx(ρ0, ρ1)

    L = M_Boson_Fermion(Hsys, tierb, tierf, Bbath, [Fbath, Fbath]; verbose=false)
    @test size(L) == (8220, 8220)
    @test L.N  == 2055
    @test L.Nb == 15
    @test L.Nf == 137
    @test nnz(L.data) == 167108
    addDissipator!(L, J)
    @test nnz(L.data) == 175330
    ados = SteadyState(L, [0.64 0; 0 0.36]; verbose=false)
    @test ados.dim == L.dim
    @test ados.Nb  == L.Nb
    @test ados.Nf  == L.Nf
    @test length(ados) == L.N
    ρ0 = ados[1, 1]
    @test getRho(ados) == ρ0
    ρ1 = [
        0.4969336345381041 - 1.190543638771177e-7im  -0.0030853446562299934 + 0.0025733625397011384im;
        -0.003085086715558734 - 0.002573561632244327im    0.5030663654832557 + 1.1905442286856098e-7im
    ]
    @test _is_Matrix_approx(ρ0, ρ1)

    ## check exceptions
    @test_throws BoundsError L[1, 8221]
    @test_throws BoundsError L[1:8221, 8220]
    @test_throws ErrorException ados[1]
    @test_throws ErrorException getADO(ados, 1)
    @test_throws ErrorException ados[L.Nb + 1, 1]
    @test_throws ErrorException ados[1, L.Nf + 1]
    @test_throws ErrorException @test_warn "Heom doesn't support matrix type : Vector{Int64}" M_Boson_Fermion([0, 0], tierb, tierf, Bbath, Fbath; verbose=false)
end

@testset "Hierarchy Dictionary" begin
    λ = 1
    W = 1
    T = 1
    μ = 1
    tier = 2
    Hsys = spzeros(ComplexF64, 2, 2)
    Q    = spzeros(ComplexF64, 2, 2)

    b1 = Boson_DrudeLorentz_Pade(Q, λ, W, T, 5)
    b2 = Boson_DrudeLorentz_Pade(Q, λ, W, T, 4)
    f1 = Fermion_Lorentz_Pade(Q, λ, μ, W, T, 3)
    f2 = Fermion_Lorentz_Pade(Q, λ, μ, W, T, 2)
    Bbath = [b1, b2]
    Fbath = [f1, f2]

    L = M_Boson_Fermion(Hsys, tier, tier, Bbath, Fbath; verbose=false)
    
    # check boson hierarchy dict.
    hDict = L.hierarchy_b
    @test L.Nb == length(hDict.idx2ado)
    for (idx, ado) in enumerate(hDict.idx2ado)
        @test hDict.ado2idx[ado] == idx
    end
    for lvl in 0:tier
        idx_list = hDict.lvl2idx[lvl]
        for idx in idx_list
            @test sum(hDict.idx2ado[idx]) == lvl
        end
    end
    for (i, b) in enumerate(Bbath)
        Ptr = hDict.bathPtr
        @test length(b) == (Ptr[i + 1] - Ptr[i])
    end

    # check fermion hierarchy dict.
    hDict = L.hierarchy_f
    @test L.Nf == length(hDict.idx2ado)
    for (idx, ado) in enumerate(hDict.idx2ado)
        @test hDict.ado2idx[ado] == idx
    end
    for lvl in 0:tier
        idx_list = hDict.lvl2idx[lvl]
        for idx in idx_list
            @test sum(hDict.idx2ado[idx]) == lvl
        end
    end
    for (i, b) in enumerate(Fbath)
        Ptr = hDict.bathPtr
        @test length(b) == (Ptr[i + 1] - Ptr[i])
    end
end
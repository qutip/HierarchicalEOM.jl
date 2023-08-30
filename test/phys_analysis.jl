# System Hamiltonian and initial state
Hsys = 0.25 * [1 0; 0 -1] + 0.5 * [0 1; 1 0]
ρ0   = [1 0; 0 0]

# Bath properties:
λ  = 0.1
W  = 0.5
kT = 0.5
N  = 2
Q = [1 0; 0 -1]  # System-bath coupling operator
bath = Boson_DrudeLorentz_Pade(Q, λ, W, kT, N)

# HEOM Liouvillian superoperator matrix
tier = 5
L = M_Boson(Hsys, tier, bath; verbose=false)

ados = SteadyState(L, ρ0; verbose=false)
ρs   = getRho(ados)
@testset "Steady state" begin
    O = [1 0.5; 0.5 1]
    @test Expect(O, ados) ≈ real(tr(O * ρs))
    @test Expect(O, ados, take_real=false) ≈ tr(O * ρs)   

    ρ1 = getRho(SteadyState(L; verbose=false))
    @test _is_Matrix_approx(ρ1, ρs)

    mat  = spzeros(ComplexF64, 2, 2)
    mat2 = spzeros(ComplexF64, 3, 3)
    bathf = Fermion_Lorentz_Pade(mat, 1, 1, 1, 1, 2)
    @test_throws ErrorException SteadyState(M_Fermion(mat, 2, bathf, :odd; verbose=false))
    @test_throws ErrorException SteadyState(M_Fermion(mat, 2, bathf, :odd; verbose=false), mat)
    @test_throws ErrorException @test_warn "The size of input matrix should be: (2, 2)." SteadyState(L, mat2)
    @test_throws ErrorException SteadyState(L, ADOs(zeros(8), 2))
    @test_throws ErrorException SteadyState(L, ADOs(ados.data, ados.N, :odd))
end

@testset "Time evolution" begin
    bath = Boson_DrudeLorentz_Pade(Q, λ, W, kT, N)
    L = M_Boson(Hsys, tier, bath; verbose=false)
    ρ0   = [1 0; 0 0]
    ρ_wrong = zeros(3, 3)

    Δt    = 10
    steps = 10
    tlist = 0:Δt:(Δt * steps)
    if isfile("evolution_p.jld2")
        rm("evolution_p.jld2")
    end
    # using the method based on propagator
    ados_list = evolution(L, ρ0, Δt, steps; verbose=false, filename="evolution_p")
    ados_wrong1 = ADOs(zeros(8), 2)
    ados_wrong2 = ADOs((ados_list[1]).data, (ados_list[1]).N, :odd)
    ρ_list_p = getRho.(ados_list)
    @test_throws ErrorException evolution(L, ρ0, Δt, steps; verbose=false, filename="evolution_p")
    @test_throws ErrorException @test_warn "The size of input matrix should be: (2, 2)." evolution(L, ρ_wrong, Δt, steps; verbose=false)
    @test_throws ErrorException evolution(L, ados_wrong1, Δt, steps)
    @test_throws ErrorException evolution(L, ados_wrong2, Δt, steps)

    if isfile("evolution_o.jld2")
        rm("evolution_o.jld2")
    end
    # using the method based on ODE solver
    ρ_list_e = getRho.(evolution(L, ρ0, tlist; verbose=false, filename="evolution_o"))
    @test_throws ErrorException evolution(L, ρ0, tlist; verbose=false, filename="evolution_o")
    @test_throws ErrorException @test_warn "The size of input matrix should be: (2, 2)." evolution(L, ρ_wrong, tlist; verbose=false)
    @test_throws ErrorException evolution(L, ados_wrong1, tlist)
    @test_throws ErrorException evolution(L, ados_wrong2, tlist)

    for i in 1:(steps + 1)
        @test _is_Matrix_approx(ρ_list_p[i], ρ_list_e[i])
    end
    @test _is_Matrix_approx(ρs, ρ_list_p[end])
    @test _is_Matrix_approx(ρs, ρ_list_e[end])

    # time-dependent Hamiltonian
    σz  = [1 0; 0 -1]
    P01 = [0 1; 0  0]

    H_sys = 0 * σz
    ρ0    = [0.5 0.5; 0.5 0.5]

    bath = Boson_DrudeLorentz_Pade(σz, 0.0005, 0.005, 0.05, 3)
    L = M_Boson(H_sys, 6, bath; verbose=false)

    function Ht(param, t)
        amplitude, delay, integral = param
        duration = integral / amplitude
        period = duration + delay
        
        t = t % period
        if t < duration
            return amplitude * [0 1; 1  0]
        else
            return 0 * zeros(2, 2)
        end
    end

    tlist = 0:10:400
    if isfile("evolution_t.jld2")
        rm("evolution_t.jld2")
    end
    fastDD_ados = evolution(L, ρ0, tlist, Ht, (0.50, 20, π/2); reltol=1e-12, abstol=1e-12, verbose=false, filename="evolution_t");
    @test_throws ErrorException evolution(L, ρ0, tlist, Ht, (0.50, 20, π/2); verbose=false, filename="evolution_t")
    @test_throws ErrorException @test_warn "The size of input matrix should be: (2, 2)." evolution(L, ρ_wrong, tlist, Ht; verbose=false)
    fastBoFiN = [
        0.4999999999999999,
        0.4972948770876402,
        0.48575366430956535,
        0.48623357850926063,
        0.4953126026707093,
        0.4956964571136953,
        0.4920357760578866,
        0.4791950252661603,
        0.4892461036798634,
        0.4928465660900118,
        0.4924652444726348,
        0.48681317030861687,
        0.4810852414658711,
        0.4895850463738865,
        0.4883063278114976,
        0.4891347504183582,
        0.4812529995344225,
        0.48397690318479186,
        0.48766721173648925,
        0.4863478506723415,
        0.4853757629933787,
        0.47619727575597826,
        0.4846116525810951,
        0.4838281128868905,
        0.4847126978506697,
        0.4809571657303493,
        0.47862209473563255,
        0.4832775448598134,
        0.4795481060217657,
        0.48232522269103417,
        0.47572829624493995,
        0.47970784844818903,
        0.48020739013048824,
        0.47927809262397914,
        0.47904939060891494,
        0.4728486106129371,
        0.47901940281020244,
        0.4755943953129941,
        0.47816314189739584,
        0.47479965067847246,
        0.47451220871416044
    ]
    fastDD = Expect(P01, fastDD_ados)
    @test typeof(fastDD) == Vector{Float64}
    for i in 1:length(tlist)
        @test fastDD[i] ≈ fastBoFiN[i] atol=1.0e-6
    end

    slowDD_ados = evolution(L, ρ0, tlist, Ht, (0.01, 20, π/2); reltol=1e-12, abstol=1e-12, verbose=false);
    slowBoFiN = [
        0.4999999999999999,
        0.4949826158957288,
        0.4808740017602084,
        0.4593383302633091,
        0.43252194834417496,
        0.402802685601788,
        0.3725210685299983,
        0.34375397372926303,
        0.3181509506463072,
        0.29684111346956915,
        0.2804084225914882,
        0.26892558062398203,
        0.26203207044640237,
        0.2590399640129812,
        0.25905156533518636,
        0.26107507001114955,
        0.2639313319708524,
        0.26368823254606255,
        0.25936367828335455,
        0.25399924039230276,
        0.2487144584444748,
        0.24360738475043517,
        0.23875449658047362,
        0.23419406109247376,
        0.2299206674661969,
        0.22588932666421335,
        0.22202631685357516,
        0.2182434870635807,
        0.21445292381955944,
        0.210579531547926,
        0.20656995199483505,
        0.20239715154725021,
        0.19806079107370503,
        0.19358407389379975,
        0.18856981068448642,
        0.18126256908840252,
        0.17215645949526584,
        0.16328689667616822,
        0.1552186906255167,
        0.14821389956195355,
        0.14240802098404504
    ]
    slowDD = Expect(P01, slowDD_ados; take_real=false)
    @test typeof(slowDD) == Vector{ComplexF64}
    for i in 1:length(tlist)
        @test slowDD[i] ≈ slowBoFiN[i] atol=1.0e-6
    end

    H_wrong1(param, t) = zeros(3, 3)
    function H_wrong2(param, t) 
        if t == 0 
            zeros(2, 2)
        else
            zeros(3, 3)
        end
    end
    @test_throws ErrorException("The dimension of `H` at t=0 is not consistent with `M.dim`.")  @test_warn "The size of input matrix should be: (2, 2)."  evolution(L, ρ0, tlist, H_wrong1; verbose=false);
    @test_throws ErrorException @test_warn "The size of input matrix should be: (2, 2)."  evolution(L, ρ0, tlist, H_wrong2; verbose=false);
    @test_throws ErrorException evolution(L, ados_wrong1, tlist, Ht)
    @test_throws ErrorException evolution(L, ados_wrong2, tlist, Ht)
end

@testset "Power spectral density" begin
    a = [0 1; 0 0]

    Hsys = a' * a

    λ  = 1e-4
    W  = 2e-1
    kT = 0.5
    N  = 5
    bath = Boson_DrudeLorentz_Matsubara((a' + a), λ, W, kT, N)

    tier = 3
    L = M_Boson(Hsys, tier, bath; verbose=false)
    L = addBosonDissipator(L, 1e-3 * a')
    L = addTerminator(L, bath)

    ados_s = SteadyState(L; verbose=false)
    ωlist = 0.9:0.01:1.1

    if isfile("PSD.txt")
        rm("PSD.txt")
    end
    psd1 = spectrum(L, ados_s, a, ωlist; verbose=false, filename="PSD")
    psd2 = [
        8.88036729e-04,
        1.06145358e-03,
        1.30081318e-03,
        1.64528197e-03,
        2.16857171e-03,
        3.02352743e-03,
        4.57224555e-03,
        7.85856353e-03,
        1.70441286e-02,
        6.49071303e-02,
        1.64934976e+02,
        6.60426108e-02,
        1.57205620e-02,
        6.74130995e-03,
        3.67130512e-03,
        2.27825666e-03,
        1.53535649e-03,
        1.09529424e-03,
        8.14605980e-04,
        6.25453434e-04,
        4.92451868e-04
    ]
    for i in 1:length(ωlist)
        @test psd1[i] ≈ psd2[i] atol=1.0e-6
    end

    mat  = spzeros(ComplexF64, 2, 2)
    mat2 = spzeros(ComplexF64, 3, 3)
    bathf = Fermion_Lorentz_Pade(mat, 1, 1, 1, 1, 2)
    @test_throws ErrorException spectrum(L, ados_s, a, ωlist; verbose=false, filename="PSD")
    @test_throws ErrorException @test_warn "The size of input matrix should be: (2, 2)." spectrum(L, ados_s, mat2, ωlist; verbose=false)
    @test_throws ErrorException spectrum(L, ADOs(zeros(8), 2), a, ωlist; verbose=false)
    @test_throws ErrorException spectrum(L, ADOs(ados_s.data, ados_s.N, :odd), a, ωlist; verbose=false)
    @test_throws ErrorException spectrum(M_Fermion(mat, 2, bathf, :odd; verbose=false), mat, mat, [0])
end

@testset "Density of states" begin
    e = -5
    U = 10
    d_up = kron(      [0 1; 0 0], [1 0; 0 1])
    d_dn = kron(-1 * [1 0; 0 -1], [0 1; 0 0])
    iden = kron(      [1 0; 0 1], [1 0; 0 1])

    H0 = e * (d_up' * d_up + d_dn' * d_dn)
    H1 = U * (d_up' * d_up * d_dn' * d_dn)
    Hsys = H0 + H1

    λ   =  1
    μ_l =  1
    μ_r = -1
    W   = 10
    kT  =  0.5
    N = 5
    fuL = Fermion_Lorentz_Pade(d_up, λ, μ_l, W, kT, N)
    fdL = Fermion_Lorentz_Pade(d_dn, λ, μ_l, W, kT, N)
    fuR = Fermion_Lorentz_Pade(d_up, λ, μ_r, W, kT, N)
    fdR = Fermion_Lorentz_Pade(d_dn, λ, μ_r, W, kT, N)

    tier = 2
    Le = M_Fermion(Hsys, tier, [fuL, fdL, fuR, fdR]; verbose=false)
    Lo = M_Fermion(Hsys, tier, [fuL, fdL, fuR, fdR], :odd; verbose=false)

    ados_s = SteadyState(Le; verbose=false)
    ωlist = -20:2:20

    if isfile("DOS.txt")
        rm("DOS.txt")
    end
    dos1 = spectrum(Lo, ados_s, d_up, ωlist; verbose=false, filename="DOS")
    dos2 = [
        0.0007920428534358747,
        0.0012795202828027256,
        0.0022148985361417936,
        0.004203651086852703,
        0.009091831276694192,
        0.024145570990254425,
        0.09111929563034224,
        0.2984323182857298,
        0.1495897397559431,
        0.12433521300531171,
        0.17217519700362036,
        0.1243352130053117,
        0.14958973975594306,
        0.29843231828572964,
        0.09111929563034214,
        0.024145570990254432,
        0.009091831276694183,
        0.004203651086852702,
        0.0022148985361417923,
        0.0012795202828027236,
        0.0007920428534358735
    ]
    for i in 1:length(ωlist)
        @test dos1[i] ≈ dos2[i] atol=1.0e-10
    end

    mat  = spzeros(ComplexF64, 2, 2)
    mat2 = spzeros(ComplexF64, 3, 3)
    bathb = Boson_DrudeLorentz_Pade(mat, 1, 1, 1, 2)
    @test_throws ErrorException spectrum(Lo, ados_s, d_up, ωlist; verbose=false, filename="DOS")
    @test_throws ErrorException @test_warn "The size of input matrix should be: (2, 2)." spectrum(Lo, ados_s, mat2, ωlist; verbose=false)
    @test_throws ErrorException spectrum(Lo, ADOs(zeros(8), 2), d_up, ωlist; verbose=false)
    @test_throws ErrorException spectrum(Lo, ADOs(ados_s.data, ados_s.N, :odd), d_up, ωlist; verbose=false)
    @test_throws ErrorException spectrum(M_Boson(mat, 2, bathb; verbose=false), mat, mat, [0])
    @test_throws ErrorException spectrum(M_Fermion(mat, 2, fuL; verbose=false), mat, mat, [0])
end

# remove all the temporary files
rm("evolution_p.jld2")
rm("evolution_o.jld2")
rm("evolution_t.jld2")
rm("PSD.txt")
rm("DOS.txt")
@testitem "Power spectrum" begin
    a = destroy(2)

    Hsys = a' * a

    λ = 1e-4
    W = 2e-1
    kT = 0.5
    N = 5
    bath = Boson_DrudeLorentz_Matsubara((a' + a), λ, W, kT, N)

    tier = 3
    L = M_Boson(Hsys, tier, bath; verbose = false)
    L = addBosonDissipator(L, 1e-3 * a')
    L = addTerminator(L, bath)

    ados_s = steadystate(L; verbose = false)
    ωlist = 0.9:0.01:1.1

    if isfile("PSD.txt")
        rm("PSD.txt")
    end
    psd1 = PowerSpectrum(L, ados_s, a, ωlist; verbose = false, filename = "PSD")
    psd2 = [
        0.0008880367286438109,
        0.0010614535771538414,
        0.0013008131768316171,
        0.00164528196789041,
        0.002168571713539654,
        0.0030235274306746876,
        0.004572245548638415,
        0.007858563532393089,
        0.017044128559126213,
        0.06490713030997987,
        164.934976224019,
        0.06604261080743092,
        0.01572056201265402,
        0.006741309950113649,
        0.003671305124632234,
        0.002278256664328282,
        0.0015353564873150376,
        0.0010952942370712003,
        0.0008146059795490796,
        0.0006254534339912124,
        0.0004924518684975365,
    ]
    for i in 1:length(ωlist)
        @test psd1[i] ≈ psd2[i]# atol = 1e-8
    end
    @test length(readlines("PSD.txt")) == length(ωlist)

    mat = Qobj(zeros(ComplexF64, 2, 2))
    mat2 = Qobj(zeros(ComplexF64, 3, 3))
    a_even = HEOMSuperOp(spre(a), EVEN, L)
    a_odd = HEOMSuperOp(spre(a), ODD, L)
    bathf = Fermion_Lorentz_Pade(mat, 1, 1, 1, 1, 2)
    @test_throws ErrorException PowerSpectrum(L, ados_s, a, ωlist; verbose = false, filename = "PSD")
    @test_throws ErrorException PowerSpectrum(L, ados_s, a_even, a_odd, ωlist; verbose = false)
    @test_throws ErrorException PowerSpectrum(L, ados_s, mat2, ωlist; verbose = false)
    @test_throws ErrorException PowerSpectrum(L, ADOs(zeros(8), 2), a, ωlist; verbose = false)
    @test_throws ErrorException PowerSpectrum(L, ADOs(zeros(32), 2), a, ωlist; verbose = false)

    # remove all the temporary files
    rm("PSD.txt")

    # two time correlation functions
    tlist = 0:0.1:1000
    corr1 = correlation_2op_1t(L, ados_s, tlist, a', a; verbose = false)
    corr2 = correlation_3op_2t(L, ados_s, [0], tlist, qeye(2), a', a; verbose = false)
    corr3 = correlation_3op_1t(L, ados_s, tlist, a', a, qeye(2), verbose = false)
    corr4 = correlation_2op_2t(L, ados_s, [0], tlist, a', a, reverse = true, verbose = false)

    @test corr1 == corr2[1, :]
    @test corr3 == corr4[1, :]
end

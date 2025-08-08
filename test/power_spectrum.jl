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
        0.0008880390384118277,
        0.0010614561417953718,
        0.0013008160598639432,
        0.0016452852599852068,
        0.0021685755505722426,
        0.0030235320297574537,
        0.00457225128899316,
        0.007858571170202312,
        0.01704413997469293,
        0.06490715292415908,
        164.93497740913972,
        0.06604258739073389,
        0.015720550418158088,
        0.0067413022487601455,
        0.003671299361034444,
        0.0022782520605674266,
        0.0015353526555966761,
        0.001095290956221469,
        0.0008146031115163797,
        0.0006254508868290593,
        0.0004924495778894673,
    ]
    for i in 1:length(ωlist)
        @test psd1[i] ≈ psd2[i] atol = 1e-8
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

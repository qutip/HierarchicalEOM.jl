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
    psd1 = PowerSpectrum(L, ados_s, a, ωlist; progress_bar = Val(true), filename = "PSD") # also test progress_bar
    psd2 = [ # result with alg = UMFPACKFactorization()
        0.0008880366931592341,
        0.0010614535356223845,
        0.001300813127188815,
        0.001645281906911075,
        0.0021685716358498205,
        0.003023527326520538,
        0.004572245398074222,
        0.007858563286641168,
        0.017044128056761452,
        0.06490712852247811,
        164.9349720231878,
        0.06604260927070985,
        0.015720561683454166,
        0.006741309825460715,
        0.0036713050661672233,
        0.002278256634144162,
        0.0015353564712408745,
        0.0010952942287540447,
        0.0008146059757785583,
        0.0006254534330017046,
        0.0004924518692553012,
    ]
    for i in 1:length(ωlist)
        @test psd1[i] ≈ psd2[i] atol = 1e-6
    end
    @test length(readlines("PSD.txt")) == length(ωlist)

    mat = Qobj(zeros(ComplexF64, 2, 2))
    mat2 = Qobj(zeros(ComplexF64, 3, 3))
    a_even = HEOMSuperOp(spre(a), EVEN, L)
    a_odd = HEOMSuperOp(spre(a), ODD, L)
    bathf = Fermion_Lorentz_Pade(mat, 1, 1, 1, 1, 2)
    @test_throws ErrorException PowerSpectrum(L, ados_s, a, ωlist; verbose = false)
    @test_throws ErrorException PowerSpectrum(L, ados_s, a, ωlist; progress_bar = Val(false), filename = "PSD")
    @test_throws ErrorException PowerSpectrum(L, ados_s, a_even, a_odd, ωlist; progress_bar = Val(false))
    @test_throws ErrorException PowerSpectrum(L, ados_s, mat2, ωlist; progress_bar = Val(false))
    @test_throws ErrorException PowerSpectrum(L, ADOs(zeros(8), 2), a, ωlist; progress_bar = Val(false))
    @test_throws ErrorException PowerSpectrum(L, ADOs(zeros(32), 2), a, ωlist; progress_bar = Val(false))

    # deprecated kwargs
    import LinearSolve # when removing this, remember to also remove LinearSolve in test/Project.toml
    @test_throws ErrorException PowerSpectrum(L, ados_s, a, ωlist; solver = LinearSolve.KrylovJL_GMRES())

    # remove all the temporary files
    rm("PSD.txt")

    # two time correlation functions
    tlist = 0:0.1:1000
    corr1 = correlation_2op_1t(L, ados_s, tlist, a', a; progress_bar = Val(false))
    corr2 = correlation_3op_2t(L, ados_s, [0], tlist, qeye(2), a', a; progress_bar = Val(false))
    corr3 = correlation_3op_1t(L, ados_s, tlist, a', a, qeye(2), progress_bar = Val(false))
    corr4 = correlation_2op_2t(L, ados_s, [0], tlist, a', a, reverse = true, progress_bar = Val(false))

    @test corr1 == corr2[1, :]
    @test corr3 == corr4[1, :]
end

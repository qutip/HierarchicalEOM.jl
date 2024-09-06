@time @testset "Power spectrum" begin
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

    ados_s = SteadyState(L; verbose = false)
    ωlist = 0.9:0.01:1.1

    if isfile("PSD.txt")
        rm("PSD.txt")
    end
    psd1 = PowerSpectrum(L, ados_s, a, ωlist; verbose = false, filename = "PSD")
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
        4.92451868e-04,
    ]
    for i in 1:length(ωlist)
        @test psd1[i] ≈ psd2[i] atol = 1.0e-6
    end

    mat = Qobj(spzeros(ComplexF64, 2, 2))
    mat2 = Qobj(spzeros(ComplexF64, 3, 3))
    a_even = HEOMSuperOp(a, EVEN, L)
    a_odd = HEOMSuperOp(a, ODD, L)
    bathf = Fermion_Lorentz_Pade(mat, 1, 1, 1, 1, 2)
    @test_throws ErrorException PowerSpectrum(L, ados_s, a, ωlist; verbose = false, filename = "PSD")
    @test_throws ErrorException PowerSpectrum(L, ados_s, a_even, ωlist; verbose = false)
    @test_throws ErrorException PowerSpectrum(L, ados_s, a_even, a_odd, ωlist; verbose = false)
    @test_throws ErrorException PowerSpectrum(L, ados_s, mat2, ωlist; verbose = false)
    @test_throws ErrorException PowerSpectrum(L, ADOs(zeros(8), 2), a, ωlist; verbose = false)
    @test_throws ErrorException PowerSpectrum(L, ADOs(zeros(32), 2), a, ωlist; verbose = false)

    # remove all the temporary files
    rm("PSD.txt")
end

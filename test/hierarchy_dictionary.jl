@time @testset "Hierarchy Dictionary" begin
    Btier = 2
    Ftier = 2
    Nb = 3
    Nf = 3
    threshold = 1e-5

    Γ = 0.0025
    Dα = 30
    Λ = 0.0025
    ωcα = 0.2
    μL = 0.5
    μR = -0.5
    kT = 0.025

    Hsys = Qobj([
        0 0 0 0
        0 0.2 0 0
        0 0 0.208 0.04
        0 0 0.04 0.408
    ])

    cop = Qobj([
        0 1 0 0
        1 0 0 0
        0 0 0 1
        0 0 1 0
    ])

    dop = Qobj([
        0 0 1 0
        0 0 0 1
        0 0 0 0
        0 0 0 0
    ])

    bbath = Boson_DrudeLorentz_Matsubara(cop, Λ, ωcα, kT, Nb)
    fbath = [Fermion_Lorentz_Pade(dop, Γ, μL, Dα, kT, Nf), Fermion_Lorentz_Pade(dop, Γ, μR, Dα, kT, Nf)]

    L = M_Boson_Fermion(Hsys, Btier, Ftier, bbath, fbath; threshold = threshold, verbose = false)
    L = addTerminator(L, bbath)

    @test size(L) == (1696, 1696)
    @test L.N == 106
    @test nnz(L.data) == 13536

    ados = SteadyState(L; verbose = false)
    @test Ic(ados, L, 1) ≈ 0.2883390125832726
    nvec_b, nvec_f = L.hierarchy.idx2nvec[1]
    @test_throws ErrorException getIndexEnsemble(nvec_f, L.hierarchy.bosonPtr)
    @test_throws ErrorException getIndexEnsemble(nvec_b, L.hierarchy.fermionPtr)
end

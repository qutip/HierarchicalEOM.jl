@testitem "HEOM superoperator" begin
    using SparseArrays

    # Take waiting time distribution as an example
    WTD_ans = [
        0.0,
        0.0015096807591478143,
        0.0019554827640519456,
        0.0019720219482432183,
        0.0018278123354472694,
        0.0016348729788834187,
        0.0014385384773227094,
        0.0012559122034857784,
        0.0010923035083891884,
        0.0009482348815344435,
        0.0008224107708326052,
        0.0007129580704171314,
        0.0006179332764127152,
        0.0005355140701828154,
        0.00046406230317384855,
        0.00040213314306247375,
        0.0003484637380485898,
        0.0003019551325507645,
        0.00026165305303984653,
        0.0002267297384879625,
        0.0001964675440128574,
    ]

    ϵ = -0.1
    Γ = 0.01
    μ = 0
    W = 1
    kT = 0.5
    N = 5
    tier = 3
    tlist = 0:50:1000
    γ_eL = γ_eR = 0.06676143449267714 # emission   rate towards Left and Right lead
    γ_aL = γ_aR = 0.0737827958503192  # absorption rate towards Left and Right lead

    d = sigmam()
    Hs = ϵ * d' * d

    ## Only local master equation
    M_me = M_S(Hs; verbose = false)
    jumpOPs = [γ_eL * d, γ_eR * d, γ_aL * d', γ_aR * d']
    M_me = addFermionDissipator(M_me, jumpOPs)

    # create quantum jump superoperator (only the emission part)
    J_me = HEOMSuperOp(γ_eR * spre(d), ODD, M_me) * HEOMSuperOp(γ_eR * spost(d'), ODD, M_me)
    @test size(M_me) == size(J_me)
    @test size(M_me, 1) == size(J_me, 1)
    @test eltype(M_me) == eltype(J_me)
    @test nnz(J_me.data) == 1
    @test J_me[4, 1] ≈ 0.0044570891355200214
    @test J_me.data * 2 ≈ (J_me - (-1) * J_me).data
    @test J_me.data * 2 ≈
        (HEOMSuperOp(√(2) * γ_eR * spre(d), ODD, M_me) * HEOMSuperOp(√(2) * γ_eR * spost(d'), ODD, M_me)).data
    @test show(devnull, MIME("text/plain"), J_me) === nothing

    # create HEOM Liouvuillian superoperator without quantum jump
    M0 = M_me - J_me
    ados_s = steadystate(M_me; verbose = false)
    ados_list = heomsolve(M0, J_me * ados_s, tlist; progress_bar = Val(false)).ados

    # calculating waiting time distribution
    WTD_me = []
    trJρ = expect(J_me, ados_s)
    WTD_me = expect(J_me, ados_list) ./ trJρ
    for i in 1:length(tlist)
        @test isapprox(WTD_me[i], WTD_ans[i]; atol = 1.0e-8)
    end
    @test isapprox(expect(J_me, ados_list[end]) / trJρ, WTD_ans[end]; atol = 1.0e-8)

    ## Left lead with HEOM method and Right lead with local master equation
    bath = Fermion_Lorentz_Pade(d, Γ, μ, W, kT, N)
    jumpOPs = [γ_eR * d, γ_aR * d']
    M_heom = M_Fermion(Hs, tier, bath; verbose = false)
    M_heom = addFermionDissipator(M_heom, jumpOPs)

    # create quantum jump superoperator (only the emission part)
    J_heom = HEOMSuperOp(γ_eR * spre(d), ODD, M_heom) * HEOMSuperOp(γ_eR * spost(d'), ODD, M_heom)

    # create HEOM Liouvuillian superoperator without quantum jump
    M1 = M_heom - J_heom
    ados_s1 = steadystate(M_heom; verbose = false)
    ados_list1 = heomsolve(M1, J_heom * ados_s1, tlist; progress_bar = Val(false)).ados

    # calculating waiting time distribution
    WTD_heom = []
    trJρ1 = expect(J_heom, ados_s1)
    WTD_heom = expect(J_heom, ados_list1) ./ trJρ1
    for i in 1:length(tlist)
        @test WTD_heom[i] ≈ WTD_ans[i] atol = 1.0e-5
    end

    ## test for error
    J_wrong1 = HEOMSuperOp(γ_eR * spre(d), ODD, M_me)
    J_wrong2 = HEOMSuperOp(γ_eR * spre(d), EVEN, M_heom)
    @test_throws ErrorException HEOMSuperOp(d, EVEN, M_heom)
    @test_throws ErrorException J_me * J_wrong2
    @test_throws ErrorException J_me + J_wrong1
    @test_throws ErrorException J_me + J_wrong2
    @test_throws ErrorException M_me + J_wrong1
    @test_throws ErrorException M_me + J_wrong2
    @test_throws ErrorException J_me - J_wrong1
    @test_throws ErrorException J_me - J_wrong2
    @test_throws ErrorException M_me - J_wrong1
    @test_throws ErrorException M_me - J_wrong2
    @test_throws ErrorException HEOMSuperOp(d, ODD, M_heom, "L")
    @test_throws ErrorException HEOMSuperOp(d, ODD, ados_s1, "L")
end

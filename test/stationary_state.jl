@time @testset "Stationary state" begin

    # System Hamiltonian and initial state
    d = sigmam()
    Hsys = d' * d
    ρ0 = Qobj([1 0; 0 0])

    # Bath properties:
    Γ = 1.0
    W = 1
    kT = 0.025851991
    μL = 1.0
    μR = -1.0
    N = 3
    baths = [Fermion_Lorentz_Pade(d, Γ, μL, W, kT, N), Fermion_Lorentz_Pade(d, Γ, μR, W, kT, N)]

    # HEOM Liouvillian superoperator matrix
    tier = 5
    L = M_Fermion(Hsys, tier, baths; verbose = false)

    ados = SteadyState(L, ρ0; verbose = false)
    ρs = getRho(ados)
    O = qeye(2) + 0.5 * sigmax()
    @test expect(O, ados) ≈ real(tr(O * ρs))
    @test expect(O, ados, take_real = false) ≈ tr(O * ρs)

    ρ1 = getRho(SteadyState(L; verbose = false))
    @test isapprox(ρ1, ρs, atol = 1e-6)

    mat = Qobj(spzeros(ComplexF64, 2, 2))
    mat2 = Qobj(spzeros(ComplexF64, 3, 3))
    bathf = Fermion_Lorentz_Pade(mat, 1, 1, 1, 1, 2)
    @test_throws ErrorException SteadyState(M_Fermion(mat, 2, bathf, ODD; verbose = false))
    @test_throws ErrorException SteadyState(M_Fermion(mat, 2, bathf, ODD; verbose = false), mat)
    @test_throws ErrorException SteadyState(L, mat2)
    @test_throws ErrorException SteadyState(L, ADOs(zeros(8), 2))
    @test_throws ErrorException SteadyState(L, ADOs(ados.data, ados.N, ODD))
    @test_throws ErrorException SteadyState(L, HEOMSuperOp(d, ODD, L) * ados)
end

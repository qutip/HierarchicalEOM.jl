# System Hamiltonian and initial state
Hsys = 0.25 * [1 0; 0 -1] + 0.5 * [0 1; 1 0]
ρ0   = [1 0; 0 0]

# Bath properties:
λ = 0.1
W = 0.5
T = 0.5
N = 2
Q = [1 0; 0 -1]  # System-bath coupling operator
bath = Boson_DrudeLorentz_Pade(Q, λ, W, T, N)

# Heom liouvillian superoperator matrix
tier = 5
L = M_Boson(Hsys, tier, bath; verbose=false)

ρs = getRho(SteadyState(L, ρ0; verbose=false))
@testset "Steady state" begin
    ρ1 = getRho(SteadyState(L; verbose=false))
    @test _is_Matrix_approx(ρ1, ρs)

    mat = spzeros(ComplexF64, 2, 2)
    bathf = Fermion_Lorentz_Pade(mat, 1, 1, 1, 1, 2)
    @test_throws ErrorException SteadyState(M_Fermion(mat, 2, bathf, :odd; verbose=false))
end

@testset "Time evolution" begin
    Δt    = 10
    steps = 10
    tlist = 0:Δt:(Δt * steps)
    ρ_list_p = getRho.(evolution(L, ρ0, Δt, steps; verbose=false))  # using the method based on propagator
    ρ_list_e = getRho.(evolution(L, ρ0, tlist; verbose=false))      # using the method based on ODE solver
    for i in 1:(steps + 1)
        @test _is_Matrix_approx(ρ_list_p[i], ρ_list_e[i])
    end
    @test _is_Matrix_approx(ρs, ρ_list_p[end])
    @test _is_Matrix_approx(ρs, ρ_list_e[end])
end

@testset "Power spectral density" begin
    a = [0 1; 0 0]

    Hsys = a' * a

    λ = 1e-4
    W = 2e-1
    T = 0.5
    N = 5
    bath = Boson_DrudeLorentz_Matsubara((a' + a), λ, W, T, N)

    tier = 3
    L = M_Boson(Hsys, tier, bath; verbose=false)
    L = addDissipator(L, 1e-3 * a')
    L = addTerminator(L, bath)

    ados_s = SteadyState(L; verbose=false)
    ωlist = 0.9:0.01:1.1
    psd1 = PSD(L, ados_s, a, ωlist; verbose=false)
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

    mat = spzeros(ComplexF64, 2, 2)
    bathf = Fermion_Lorentz_Pade(mat, 1, 1, 1, 1, 2)
    @test_throws ErrorException PSD(M_Fermion(mat, 2, bathf, :odd; verbose=false), mat, mat, [0])
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

    λ = 1
    μ_l =  1
    μ_r = -1
    W = 10
    T = 0.5
    N = 5
    fuL = Fermion_Lorentz_Pade(d_up, λ, μ_l, W, T, N)
    fdL = Fermion_Lorentz_Pade(d_dn, λ, μ_l, W, T, N)
    fuR = Fermion_Lorentz_Pade(d_up, λ, μ_r, W, T, N)
    fdR = Fermion_Lorentz_Pade(d_dn, λ, μ_r, W, T, N)

    tier = 2
    Le = M_Fermion(Hsys, tier, [fuL, fdL, fuR, fdR]; verbose=false)
    Lo = M_Fermion(Hsys, tier, [fuL, fdL, fuR, fdR], :odd; verbose=false)

    ados_s = SteadyState(Le; verbose=false)
    ωlist = -20:2:20
    dos1 = DOS(Lo, ados_s, d_up, ωlist; verbose=false)
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

    mat = spzeros(ComplexF64, 2, 2)
    bathb = Boson_DrudeLorentz_Pade(mat, 1, 1, 1, 2)
    @test_throws ErrorException DOS(M_Boson(mat, 2, bathb; verbose=false), mat, mat, [0])
    @test_throws ErrorException DOS(M_Fermion(mat, 2, fuL; verbose=false), mat, mat, [0])
end
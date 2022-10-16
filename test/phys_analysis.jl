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
    Me = M_Fermion(Hsys, tier, [fuL, fdL, fuR, fdR]; verbose=false)
    Mo = M_Fermion(Hsys, tier, [fuL, fdL, fuR, fdR], :odd; verbose=false)

    ados_s = SteadyState(Me; verbose=false)
    ωlist = 0:2:20
    dos1 = DOS(Mo, ados_s, d_up, ωlist; verbose=false)
    dos2 = [
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
end
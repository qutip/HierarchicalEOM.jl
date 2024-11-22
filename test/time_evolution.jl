@time @testset "Time evolution" begin

    # System Hamiltonian and initial state
    Hsys = 0.25 * sigmaz() + 0.5 * sigmax()
    ψ0 = basis(2, 0)

    # Bath properties:
    λ = 0.1
    W = 0.5
    kT = 0.5
    N = 2
    tier = 5
    Q = sigmaz()  # System-bath coupling operator
    e_ops = [rand_dm(2), rand_dm(2)]

    bath = Boson_DrudeLorentz_Pade(Q, λ, W, kT, N)

    L = M_Boson(Hsys, tier, bath; verbose = false)
    ρs = getRho(steadystate(L; verbose = false))
    ρ_wrong = Qobj(zeros(3, 3))

    Δt = 10
    steps = 10
    tlist = 0:Δt:(Δt*steps)
    if isfile("evolution_p.jld2")
        rm("evolution_p.jld2")
    end
    # using the method based on propagator
    ados_list = HEOMsolve(L, ψ0, Δt, steps; verbose = false, filename = "evolution_p").ados
    sol_p = HEOMsolve(L, ψ0, Δt, steps; e_ops = e_ops, verbose = false)
    expvals_p = sol_p.expect
    ados_wrong1 = ADOs(zeros(8), 2)
    ados_wrong2 = ADOs(zeros(32), 2)
    ados_wrong3 = ADOs((ados_list[1]).data, (ados_list[1]).N, ODD)
    ados_wrong4 = HEOMSuperOp(spre(Q), ODD, ados_list[end]) * ados_list[end]
    ρ_list_p = getRho.(ados_list)
    @test show(devnull, MIME("text/plain"), sol_p) === nothing
    @test length(sol_p.ados) == 1
    @test_throws ErrorException evolution(L, ψ0, Δt, steps; verbose = false)
    @test_throws ErrorException HEOMsolve(L, ψ0, Δt, steps; verbose = false, filename = "evolution_p")
    @test_throws ErrorException HEOMsolve(L, ρ_wrong, Δt, steps; verbose = false)
    @test_throws ErrorException HEOMsolve(L, ados_wrong1, Δt, steps; verbose = false)
    @test_throws ErrorException HEOMsolve(L, ados_wrong2, Δt, steps; verbose = false)
    @test_throws ErrorException HEOMsolve(L, ados_wrong3, Δt, steps; verbose = false)
    @test_throws ErrorException HEOMsolve(L, ados_wrong4, Δt, steps; verbose = false)

    if isfile("evolution_o.jld2")
        rm("evolution_o.jld2")
    end
    # using the method based on ODE solver
    sol_e = HEOMsolve(L, ψ0, tlist; e_ops = e_ops, saveat = tlist, verbose = false, filename = "evolution_o")
    ρ_list_e = getRho.(sol_e.ados)
    expvals_e = sol_e.expect
    @test show(devnull, MIME("text/plain"), sol_e) === nothing
    @test_throws ErrorException evolution(L, ψ0, tlist; verbose = false)
    @test_throws ErrorException HEOMsolve(L, ψ0, tlist; verbose = false, filename = "evolution_o")
    @test_throws ErrorException HEOMsolve(L, ρ_wrong, tlist; verbose = false)
    @test_throws ErrorException HEOMsolve(L, ados_wrong1, tlist; verbose = false)
    @test_throws ErrorException HEOMsolve(L, ados_wrong2, tlist; verbose = false)
    @test_throws ErrorException HEOMsolve(L, ados_wrong3, tlist; verbose = false)
    @test_throws ErrorException HEOMsolve(L, ados_wrong4, tlist; verbose = false)

    @test all(expvals_p .≈ expvals_e)
    @test all([ρ_list_p[i] ≈ ρ_list_e[i] for i in 1:(steps+1)])
    @test isapprox(ρs, ρ_list_p[end]; atol = 1e-4)
    @test isapprox(ρs, ρ_list_e[end]; atol = 1e-4)
    @test isapprox(ρs, getRho(sol_p.ados[1]); atol = 1e-4)
    jldopen("evolution_p.jld2", "r") do file
        ados_list = file["ados"]
        @test typeof(ados_list) == Vector{ADOs}
        @test length(ados_list) == steps + 1
        @test all([ρ_list_p[i] ≈ getRho(ados_list[i]) for i in 1:(steps+1)])
    end
    jldopen("evolution_o.jld2", "r") do file
        ados_list = file["ados"]
        @test typeof(ados_list) == Vector{ADOs}
        @test length(ados_list) == steps + 1
        @test all([ρ_list_e[i] ≈ getRho(ados_list[i]) for i in 1:(steps+1)])
    end

    # time-dependent Hamiltonian
    σz = sigmaz()
    P01 = basis(2, 0) * basis(2, 1)'

    H_sys = 0 * σz
    ψ0 = (basis(2, 0) + basis(2, 1)) / √2

    bath = Boson_DrudeLorentz_Pade(σz, 0.0005, 0.005, 0.05, 3)
    L = M_Boson(H_sys, 6, bath; verbose = false)

    function Ht(t, p)
        duration = p.integral / p.amplitude
        period = duration + p.delay

        t = t % period
        if t < duration
            return p.amplitude * sigmax()
        else
            return Qobj([0 0; 0 0])
        end
    end

    tlist = 0:10:400
    if isfile("evolution_t.jld2")
        rm("evolution_t.jld2")
    end
    p_fast = (amplitude = 0.5, delay = 20, integral = π / 2)
    fastDD_sol = HEOMsolve(
        L,
        ψ0,
        tlist;
        H_t = Ht,
        params = p_fast,
        e_ops = [P01],
        saveat = tlist,
        reltol = 1e-12,
        abstol = 1e-12,
        verbose = false,
        filename = "evolution_t",
    )
    fastDD_ados = fastDD_sol.ados
    fastDD1 = real.(fastDD_sol.expect[1, :])
    fastDD2 = expect(P01, fastDD_ados)
    jldopen("evolution_t.jld2", "r") do file
        ados_list = file["ados"]
        @test typeof(ados_list) == Vector{ADOs}
        @test length(ados_list) == length(tlist)
    end
    @test_throws ErrorException HEOMsolve(
        L,
        ψ0,
        tlist;
        H_t = Ht,
        params = p_fast,
        verbose = false,
        filename = "evolution_t",
    )
    @test_throws ErrorException HEOMsolve(L, ρ_wrong, tlist; H_t = Ht, verbose = false)
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
        0.47451220871416044,
    ]
    @test show(devnull, MIME("text/plain"), fastDD_sol) === nothing
    @test length(fastDD_sol.ados) == length(tlist)
    @test size(fastDD_sol.expect) == (1, length(tlist))
    @test typeof(fastDD1) == typeof(fastDD2) == Vector{Float64}
    @test all(isapprox.(fastDD1, fastBoFiN; atol = 1.0e-6))
    @test all(isapprox.(fastDD2, fastBoFiN; atol = 1.0e-6))

    p_slow = (amplitude = 0.01, delay = 20, integral = π / 2)
    slowDD_sol = HEOMsolve(
        L,
        ψ0,
        tlist;
        H_t = Ht,
        params = p_slow,
        e_ops = [P01],
        saveat = tlist,
        reltol = 1e-12,
        abstol = 1e-12,
        verbose = false,
    )
    slowDD_ados = slowDD_sol.ados
    slowDD1 = slowDD_sol.expect[1, :]
    slowDD2 = expect(P01, slowDD_ados; take_real = false)
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
        0.14240802098404504,
    ]
    @test show(devnull, MIME("text/plain"), slowDD_sol) === nothing
    @test length(slowDD_sol.ados) == length(tlist)
    @test size(slowDD_sol.expect) == (1, length(tlist))
    @test typeof(slowDD1) == typeof(slowDD2) == Vector{ComplexF64}
    @test all(isapprox.(slowDD1, slowBoFiN; atol = 1.0e-6))
    @test all(isapprox.(slowDD2, slowBoFiN; atol = 1.0e-6))

    H_wrong1(t, p) = Qobj(zeros(3, 3))
    H_wrong2(t, p) = t == 0 ? Qobj(zeros(2, 2)) : Qobj(zeros(3, 3))
    ados_wrong1 = ADOs(zeros(8), 2)
    ados_wrong2 = ADOs(zeros(32), 2)
    ados_wrong3 = ADOs((slowDD_ados[1]).data, (slowDD_ados[1]).N, ODD)
    @test_throws ErrorException HEOMsolve(L, ψ0, tlist; H_t = H_wrong1, verbose = false)
    @test_throws ErrorException HEOMsolve(L, ψ0, tlist; H_t = H_wrong2, verbose = false)
    @test_throws ErrorException HEOMsolve(L, ados_wrong1, tlist; H_t = Ht, verbose = false)
    @test_throws ErrorException HEOMsolve(L, ados_wrong2, tlist; H_t = Ht, verbose = false)
    @test_throws ErrorException HEOMsolve(L, ados_wrong3, tlist; H_t = Ht, verbose = false)

    # remove all the temporary files
    rm("evolution_p.jld2")
    rm("evolution_o.jld2")
    rm("evolution_t.jld2")
end

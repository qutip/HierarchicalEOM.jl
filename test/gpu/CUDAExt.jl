CUDA.@time @testset "CUDA Extension" begin

    # re-define the bath (make the matrix smaller)
    λ = 0.01
    W = 0.5
    kT = 0.5
    μ = 0
    N = 3
    tier = 3

    # System Hamiltonian
    Hsys = Qobj(zeros(ComplexF64, 2, 2))

    # system-bath coupling operator
    Qb = sigmax()
    Qf = sigmam()

    E = Qobj(rand(ComplexF64, 2, 2))
    e_ops = [E]

    # initial state
    ψ0 = basis(2, 1)

    Bbath = Boson_DrudeLorentz_Pade(Qb, λ, W, kT, N)
    Fbath = Fermion_Lorentz_Pade(Qf, λ, μ, W, kT, N)

    # Solving time Evolution
    ## Schrodinger HEOMLS
    L_cpu = M_S(Hsys; verbose = false)
    L_gpu = cu(L_cpu)
    sol_cpu = HEOMsolve(L_cpu, ψ0, [0, 10]; e_ops = e_ops, verbose = false)
    sol_gpu = HEOMsolve(L_gpu, ψ0, [0, 10]; e_ops = e_ops, verbose = false)
    @test all(isapprox.(sol_cpu.expect[1, :], sol_gpu.expect[1, :], atol = 1e-4))
    @test isapprox(getRho(sol_cpu.ados[end]), getRho(sol_gpu.ados[end]), atol = 1e-4)

    ## Boson HEOMLS
    L_cpu = M_Boson(Hsys, tier, Bbath; verbose = false)
    L_gpu = cu(L_cpu)
    sol_cpu = HEOMsolve(L_cpu, ψ0, [0, 10]; e_ops = e_ops, verbose = false)
    sol_gpu = HEOMsolve(L_gpu, ψ0, [0, 10]; e_ops = e_ops, verbose = false)
    @test all(isapprox.(sol_cpu.expect[1, :], sol_gpu.expect[1, :], atol = 1e-4))
    @test isapprox(getRho(sol_cpu.ados[end]), getRho(sol_gpu.ados[end]), atol = 1e-4)

    ## Fermion HEOMLS
    L_cpu = M_Fermion(Hsys, tier, Fbath; verbose = false)
    L_gpu = cu(L_cpu)
    sol_cpu = HEOMsolve(L_cpu, ψ0, [0, 10]; e_ops = e_ops, verbose = false)
    sol_gpu = HEOMsolve(L_gpu, ψ0, [0, 10]; e_ops = e_ops, verbose = false)
    @test all(isapprox.(sol_cpu.expect[1, :], sol_gpu.expect[1, :], atol = 1e-4))
    @test isapprox(getRho(sol_cpu.ados[end]), getRho(sol_gpu.ados[end]), atol = 1e-4)

    ## Boson Fermion HEOMLS
    L_cpu = M_Boson_Fermion(Hsys, tier, tier, Bbath, Fbath; verbose = false)
    L_gpu = cu(L_cpu)
    tlist = 0:1:10
    sol_cpu = HEOMsolve(L_cpu, ψ0, tlist; e_ops = e_ops, saveat = tlist, verbose = false)
    sol_gpu = HEOMsolve(L_gpu, ψ0, tlist; e_ops = e_ops, saveat = tlist, verbose = false)
    @test all(isapprox.(sol_cpu.expect[1, :], sol_gpu.expect[1, :], atol = 1e-4))
    for i in 1:length(tlist)
        @test isapprox(getRho(sol_cpu.ados[i]), getRho(sol_gpu.ados[i]), atol = 1e-4)
    end

    # SIAM
    ϵ = -5
    U = 10
    σm = sigmam() ## σ-
    σz = sigmaz() ## σz
    II = qeye(2)  ## identity matrix
    d_up = tensor(σm, II)
    d_dn = tensor(-1 * σz, σm)
    ψ0 = tensor(basis(2, 0), basis(2, 0))
    Hsys = ϵ * (d_up' * d_up + d_dn' * d_dn) + U * (d_up' * d_up * d_dn' * d_dn)
    Γ = 2
    μ = 0
    W = 10
    kT = 0.5
    N = 5
    tier = 3
    bath_up = Fermion_Lorentz_Pade(d_up, Γ, μ, W, kT, N)
    bath_dn = Fermion_Lorentz_Pade(d_dn, Γ, μ, W, kT, N)
    bath_list = [bath_up, bath_dn]

    ## solve stationary state
    L_even_cpu = M_Fermion(Hsys, tier, bath_list; verbose = false)
    L_even_gpu = cu(L_even_cpu)
    ados_cpu = steadystate(L_even_cpu; verbose = false)
    ados_gpu = steadystate(L_even_gpu, ψ0, 10; verbose = false)
    @test all(isapprox.(ados_cpu.data, ados_gpu.data; atol = 1e-6))

    ## solve density of states
    ωlist = -5:0.5:5
    L_odd_cpu = M_Fermion(Hsys, tier, bath_list, ODD; verbose = false)
    L_odd_gpu = cu(L_odd_cpu)
    dos_cpu = DensityOfStates(L_odd_cpu, ados_cpu, d_up, ωlist; verbose = false)
    dos_gpu = DensityOfStates(
        L_odd_gpu,
        ados_cpu,
        d_up,
        ωlist;
        solver = KrylovJL_BICGSTAB(rtol = 1.0f-10, atol = 1.0f-12),
        verbose = false,
    )
    for (i, ω) in enumerate(ωlist)
        @test dos_cpu[i] ≈ dos_gpu[i] atol = 1e-6
    end
end

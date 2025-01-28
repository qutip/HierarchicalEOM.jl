@time @testset "Input-Output" begin
    tier = 12
    Δ = 2 * π
    Γ = 0.1 * Δ
    λ = 0.1 * Δ
    ω0 = 0.2 * Δ
    Z = sigmaz()
    X = sigmax()
    Hsys = 0.5 * Δ * Z
    ρ0 = ket2dm(basis(2, 0))
    bath = BosonBath(
        X,
        [0.5 * λ^2, 0.5 * λ^2],
        [-1.0im * ω0 + Γ, 1.0im * ω0 + Γ],
        [0.5im * λ^2, -0.5im * λ^2],
        [-1.0im * ω0 + Γ, 1.0im * ω0 + Γ],
    )

    # compare result with mesolve
    a = qeye(2) ⊗ destroy(tier)
    H = Hsys ⊗ qeye(tier) + λ * tensor(X, qeye(tier)) * (a + a') + ω0 * a' * a
    Ntime = 100
    tlist = LinRange(0, 20 / Δ, Ntime)
    sol_me0 = mesolve(
        H,
        ρ0 ⊗ ket2dm(basis(tier, 0)),
        tlist,
        [sqrt(Γ * 2) * a],
        e_ops = [Z ⊗ qeye(tier), a' * a],
        progress_bar = false,
    )
    sol_me1 = mesolve(
        H,
        ρ0 ⊗ ket2dm(basis(tier, 1)),
        tlist,
        [sqrt(Γ * 2) * a],
        e_ops = [Z ⊗ qeye(tier), a' * a],
        progress_bar = false,
    )

    # dynamical field
    input1(p, t) = λ * exp(-1.0im * ω0 * t - Γ * t)
    input2(p, t) = λ * exp(1.0im * ω0 * t - Γ * t)
    bath_input = BosonDynamicalField(X, η_in = [input1, input2])

    bath_output_1R = BosonDynamicalField(X, η_out_R = [λ], γ_out_R = [-1.0im * ω0 + Γ])
    bath_output_2L = BosonDynamicalField(X, η_out_L = [λ], γ_out_L = [1.0im * ω0 + Γ])

    output_1R(p, t) = λ * exp(1.0im * ω0 * (p.tout - t) - Γ * (p.tout - t))
    bath_output_1R_fn = BosonDynamicalField(X, η_out_fn_R = output_1R)

    output_2L(p, t) = λ * exp(-1.0im * ω0 * (p.tout - t) - Γ * (p.tout - t))
    bath_output_2L_fn = BosonDynamicalField(X, η_out_fn_L = output_2L)

    # only input
    bath_list_in = [bath, bath_input]
    M_in = M_Boson(Hsys, tier, bath_list_in, verbose = false)
    HDict_in = M_in.hierarchy
    Z_super = HEOMSuperOp(spre(Z), EVEN, M_in)
    e_ops_in =
        [Tr_ADO(M_in, 1) * Z_super, (Tr_ADO(M_in, 1) - Tr_ADO(M_in, HDict_in.nvec2idx[Nvec([0, 0, 1, 1])])) * Z_super]
    sol_heom_in = HEOMsolve(M_in, ρ0, tlist, e_ops = e_ops_in, verbose = false)
    @test length(sol_heom_in.ados) == 1
    @test size(sol_heom_in.expect) == (length(e_ops_in), Ntime)
    @test all(isapprox.(sol_heom_in.expect[1, :], sol_me0.expect[1, :], atol = 1e-6))
    @test all(isapprox.(sol_heom_in.expect[2, :], sol_me1.expect[1, :], atol = 1e-6))

    # only output
    bath_list_out = [bath, bath_output_1R, bath_output_2L, bath_input]
    M_out = M_Boson(Hsys, tier, bath_list_out, verbose = false)
    HDict_out = M_out.hierarchy
    e_ops_out = [-Tr_ADO(M_out, HDict_out.nvec2idx[Nvec([0, 0, 1, 1, 0, 0])])]
    sol_heom_out = HEOMsolve(M_out, ρ0, tlist, e_ops = e_ops_out, verbose = false)
    @test length(sol_heom_out.ados) == 1
    @test size(sol_heom_out.expect) == (length(e_ops_out), Ntime)
    @test all(isapprox.(sol_heom_out.expect[1, :], sol_me0.expect[2, :], atol = 1e-6))

    # output conditional on input
    e_ops_out_in = [
        Tr_ADO(M_out, 1),
        Tr_ADO(M_out, HDict_out.nvec2idx[Nvec([0, 0, 0, 1, 1, 0])]),
        Tr_ADO(M_out, HDict_out.nvec2idx[Nvec([0, 0, 1, 0, 0, 1])]),
        Tr_ADO(M_out, HDict_out.nvec2idx[Nvec([0, 0, 1, 1, 0, 0])]),
        Tr_ADO(M_out, HDict_out.nvec2idx[Nvec([0, 0, 1, 1, 1, 1])]),
    ]
    sol_heom_out_in = HEOMsolve(M_out, ρ0, tlist, e_ops = e_ops_out_in, verbose = false)
    @test length(sol_heom_out.ados) == 1
    @test size(sol_heom_out_in.expect) == (length(e_ops_out_in), Ntime)
    for (i, t) in enumerate(tlist)
        @test isapprox(
            sol_me1.expect[2, i],
            (
                exp(-2 * Γ * t) * sol_heom_out_in.expect[1, i] -
                exp(1im * ω0 * t - Γ * t) * sol_heom_out_in.expect[2, i] -
                exp(-1im * ω0 * t - Γ * t) * sol_heom_out_in.expect[3, i] - sol_heom_out_in.expect[4, i] +
                sol_heom_out_in.expect[5, i]
            ),
            atol = 1e-6,
        )
    end

    # time-dependent function output conditional on input
    bath_list_out_fn = [bath, bath_output_1R_fn, bath_output_2L_fn, bath_input]
    M_out_fn = M_Boson(Hsys, tier, bath_list_out_fn, verbose = false)
    HDict_out_fn = M_out_fn.hierarchy
    e_ops_out_fn = [
        Tr_ADO(M_out_fn, 1),
        Tr_ADO(M_out_fn, HDict_out_fn.nvec2idx[Nvec([0, 0, 0, 1, 1, 0])]),
        Tr_ADO(M_out_fn, HDict_out_fn.nvec2idx[Nvec([0, 0, 1, 0, 0, 1])]),
        Tr_ADO(M_out_fn, HDict_out_fn.nvec2idx[Nvec([0, 0, 1, 1, 0, 0])]),
        Tr_ADO(M_out_fn, HDict_out_fn.nvec2idx[Nvec([0, 0, 1, 1, 1, 1])]),
    ]
    for (i, tout) in enumerate(tlist)
        (i == 1) && continue ## TODO: fix this
        p = (tout = tout,)
        sol_heom_out_fn = HEOMsolve(M_out_fn, ρ0, [0, tout], params = p, e_ops = e_ops_out_fn, verbose = false)
        @test length(sol_heom_out_fn.ados) == 1
        @test size(sol_heom_out_fn.expect) == (length(e_ops_out_fn), 2)
        @test isapprox(
            sol_me1.expect[2, i],
            (
                exp(-2 * Γ * tout) * sol_heom_out_fn.expect[1, end] -
                exp(1im * ω0 * tout - Γ * tout) * sol_heom_out_fn.expect[2, end] -
                exp(-1im * ω0 * tout - Γ * tout) * sol_heom_out_fn.expect[3, end] - sol_heom_out_fn.expect[4, end] +
                sol_heom_out_fn.expect[5, end]
            ),
            atol = 1e-6,
        )
    end
end

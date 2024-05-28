@time @testset "Bath and Exponent" begin

    # prepare coupling operator and coefficients of exponential-exponential-expansion terms
    η0 = [1]
    γ0 = [2]
    η1 = [1, 3, 5, 7, 9]
    η2 = [2, 4, 6, 8, 10]
    γ1 = [0.1, 0.3, 0.5, 0.3, 0.7]
    γ2 = [0.1, 0.2, 0.5, 0.6, 0.9]
    γ3 = [0.1, 0.2 + 0.3im, 0.4 - 0.5im, 0.6 + 0.7im, 0.8 - 0.9im]
    γ4 = [0.1, 0.2 - 0.3im, 0.4 + 0.5im, 0.6 - 0.7im, 0.8 + 0.9im]
    op = [0 0; 0 0]

    ################################################
    # Boson bath
    b = BosonBath(op, η1, γ1, combine = false)
    @test length(b) == 5

    ## check for combine
    b = BosonBath(op, η1, γ1)
    @test length(b) == 4

    ## check for η and γ list, and coupling operator
    η = []
    γ = []
    for e in b
        push!(η, e.η)
        push!(γ, e.γ)
        @test e.op == [0 0; 0 0]
    end
    @test η == [1.0 + 0.0im, 10.0 + 0.0im, 5.0 + 0.0im, 9.0 + 0.0im]
    @test γ == [0.1 + 0.0im, 0.3 + 0.0im, 0.5 + 0.0im, 0.7 + 0.0im]

    ## check for real and image sperarate case
    bs = BosonBath(op, η1, γ1, η2, γ2, combine = false)
    @test length(bs) == 10

    ## check for combine
    b = BosonBath(op, η1, γ1, η2, γ2)
    @test length(b) == 7
    @test C(b, [0.183183])[1] ≈ C(bs, [0.183183])[1]

    ## check for η and γ list, and coupling operator
    η = []
    γ = []
    for e in b
        push!(η, e.η)
        push!(γ, e.γ)
        @test e.op == [0 0; 0 0]
        @test (e.types == "bR") || (e.types == "bI") || (e.types == "bRI")
    end
    @test η == [10.0 + 0.0im, 9.0 + 0.0im, 4.0 + 0.0im, 8.0 + 0.0im, 10.0 + 0.0im, 1.0 + 2.0im, 5.0 + 6.0im]
    @test γ == [0.3 + 0.0im, 0.7 + 0.0im, 0.2 + 0.0im, 0.6 + 0.0im, 0.9 + 0.0im, 0.1 + 0.0im, 0.5 + 0.0im]
    @test show(devnull, MIME("text/plain"), b) == nothing

    ## check for exponents
    @test show(devnull, MIME("text/plain"), b[1]) == nothing
    @test show(devnull, MIME("text/plain"), b[:]) == nothing
    @test show(devnull, MIME("text/plain"), b[1:end]) == nothing

    ## check exceptions
    @test_throws ErrorException BosonBath(op, [0], [0, 0])
    @test_throws ErrorException BosonBath(op, [0, 0], [0, 0], [0], [0, 0])
    @test_throws ErrorException BosonBath(op, [0, 0], [0], [0, 0], [0, 0])
    @test_throws ErrorException BosonBath([0, 0], [0, 0], [0, 0], [0, 0], [0, 0])
    @test_warn "The system-bosonic-bath coupling operator \"op\" should be Hermitian operator" BosonBath(
        [0 1; 0 0],
        [0, 0],
        [0, 0],
        [0, 0],
        [0, 0],
    )
    ################################################

    ################################################
    # Boson bath (RWA)
    b = BosonBathRWA(op, η1, γ3, η2, γ4)
    @test length(bs) == 10
    cp, cm = C(b, [0.183183])
    @test cp[1] ≈ 22.384506765987076 + 0.7399082821797519im
    @test cm[1] ≈ 26.994911851482776 - 0.799138487523946im

    ## check for η and γ list, and coupling operator
    η = []
    γ = []
    for e in b
        push!(η, e.η)
        push!(γ, e.γ)
        @test e.op == [0 0; 0 0]
        @test (e.types == "bA") || (e.types == "bE")
    end
    @test η == [1, 3, 5, 7, 9, 2, 4, 6, 8, 10]
    @test γ == [
        0.1,
        0.2 + 0.3im,
        0.4 - 0.5im,
        0.6 + 0.7im,
        0.8 - 0.9im,
        0.1,
        0.2 - 0.3im,
        0.4 + 0.5im,
        0.6 - 0.7im,
        0.8 + 0.9im,
    ]
    @test show(devnull, MIME("text/plain"), b) == nothing

    ## check for exponents
    @test show(devnull, MIME("text/plain"), b[1]) == nothing
    @test show(devnull, MIME("text/plain"), b[:]) == nothing
    @test show(devnull, MIME("text/plain"), b[1:end]) == nothing

    ## check exceptions
    @test_throws ErrorException BosonBathRWA(op, [0], [0, 0], [0, 0], [0, 0])
    @test_throws ErrorException BosonBathRWA(op, [0, 0], [0], [0, 0], [0, 0])
    @test_throws ErrorException BosonBathRWA(op, [0, 0], [0, 0], [0], [0, 0])
    @test_throws ErrorException BosonBathRWA(op, [0, 0], [0, 0], [0, 0], [0])
    @test_throws ErrorException BosonBathRWA([0, 0], [0, 0], [0, 0], [0, 0], [0, 0])
    @test_warn "The elements in \'γ_absorb\' should be complex conjugate of the corresponding elements in \'γ_emit\'." BosonBathRWA(
        op,
        [0, 0],
        [0, 1im],
        [0, 0],
        [0, 1im],
    )
    ################################################

    ################################################
    # Fermion bath
    b = FermionBath(op, η1, γ3, η2, γ4)
    @test length(b) == 10
    cp, cm = C(b, [0.183183])
    @test cp[1] ≈ 22.384506765987076 + 0.7399082821797519im
    @test cm[1] ≈ 26.994911851482776 - 0.799138487523946im
    for e in b
        @test e.op == [0 0; 0 0]
        @test (e.types == "fA") || (e.types == "fE")
    end
    @test show(devnull, MIME("text/plain"), b) == nothing

    ## check for exponents
    @test show(devnull, MIME("text/plain"), b[1]) == nothing
    @test show(devnull, MIME("text/plain"), b[:]) == nothing
    @test show(devnull, MIME("text/plain"), b[1:end]) == nothing

    ## check exceptions
    @test_throws ErrorException FermionBath(op, [0], [0, 0], [0, 0], [0, 0])
    @test_throws ErrorException FermionBath(op, [0, 0], [0], [0, 0], [0, 0])
    @test_throws ErrorException FermionBath(op, [0, 0], [0, 0], [0], [0, 0])
    @test_throws ErrorException FermionBath(op, [0, 0], [0, 0], [0, 0], [0])
    @test_throws ErrorException FermionBath([0, 0], [0, 0], [0, 0], [0, 0], [0, 0])
    @test_warn "The elements in \'γ_absorb\' should be complex conjugate of the corresponding elements in \'γ_emit\'." FermionBath(
        op,
        [0, 0],
        [0, 1im],
        [0, 0],
        [0, 1im],
    )
    ################################################

    ################################################
    # Exponent
    @test C(BosonBath(op, η0, γ0), [0.123])[1] ≈ η0[1] * exp(-γ0[1] * 0.123)

    ## check exceptions
    @test_throws ErrorException b[11]
    @test_throws ErrorException b[1:11]
    @test_throws ErrorException b[0:10]
    ################################################
end

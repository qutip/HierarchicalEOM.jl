@testitem "Bath and Exponent" begin

    # prepare coupling operator and coefficients of exponential-exponential-expansion terms
    η0 = [1]
    γ0 = [2]
    η1 = [1, 3, 5, 7, 9]
    η2 = [2, 4, 6, 8, 10]
    γ1 = [0.1, 0.3, 0.5, 0.3, 0.7]
    γ2 = [0.1, 0.2, 0.5, 0.6, 0.9]
    γ3 = [0.1, 0.2 + 0.3im, 0.4 - 0.5im, 0.6 + 0.7im, 0.8 - 0.9im]
    γ4 = [0.1, 0.2 - 0.3im, 0.4 + 0.5im, 0.6 - 0.7im, 0.8 + 0.9im]
    op = Qobj([0 0; 0 0])

    ################################################
    # Boson bath
    b = BosonBath(op, η1, γ1, combine = false)
    @test length(b) == 5

    ## check for combine
    b = BosonBath(op, η1, γ1)
    types = ["bRI", "bRI", "bRI", "bRI"]
    @test length(b) == 4

    ## check for η and γ list, and coupling operator
    η = Vector{ComplexF64}(undef, length(b))
    γ = Vector{ComplexF64}(undef, length(b))
    for (i, e) in enumerate(b)
        η[i] = e.η
        γ[i] = e.γ
        @test e.op == op
        @test e.types == types[i]
    end
    @test η == [1.0 + 0.0im, 10.0 + 0.0im, 5.0 + 0.0im, 9.0 + 0.0im]
    @test γ == [0.1 + 0.0im, 0.3 + 0.0im, 0.5 + 0.0im, 0.7 + 0.0im]

    ## check for real and image sperarate case
    bs = BosonBath(op, η1, γ1, η2, γ2, combine = false)
    @test length(bs) == 10

    ## check for combine
    b = BosonBath(op, η1, γ1, η2, γ2)
    types = ["bRI", "bRI", "bR", "bR", "bI", "bI", "bI"]
    @test length(b) == 7
    @test correlation_function(b, 0.183183)[1] ≈ correlation_function(bs, [0.183183])[1]
    @test_throws ErrorException C(bs, [0.183183])

    ## check for η and γ list, and coupling operator
    η = Vector{ComplexF64}(undef, length(b))
    γ = Vector{ComplexF64}(undef, length(b))
    for (i, e) in enumerate(b)
        η[i] = e.η
        γ[i] = e.γ
        @test e.op == op
        @test e.types == types[i]
    end
    @test η == [1.0 + 2.0im, 5.0 + 6.0im, 10.0 + 0.0im, 9.0 + 0.0im, 4.0 + 0.0im, 8.0 + 0.0im, 10.0 + 0.0im]
    @test γ == [0.1 + 0.0im, 0.5 + 0.0im, 0.3 + 0.0im, 0.7 + 0.0im, 0.2 + 0.0im, 0.6 + 0.0im, 0.9 + 0.0im]
    @test show(devnull, MIME("text/plain"), b) === nothing

    ## check for exponents
    @test show(devnull, MIME("text/plain"), b[1]) === nothing
    @test show(devnull, MIME("text/plain"), b[:]) === nothing
    @test show(devnull, MIME("text/plain"), b[1:end]) === nothing

    ## check exceptions
    @test_throws ErrorException BosonBath(op, [0], [0, 0])
    @test_throws ErrorException BosonBath(op, [0, 0], [0, 0], [0], [0, 0])
    @test_throws ErrorException BosonBath(op, [0, 0], [0], [0, 0], [0, 0])
    @test_throws ErrorException BosonBath(Qobj([0, 0]), [0, 0], [0, 0], [0, 0], [0, 0])
    @test_warn "The system-bosonic-bath coupling operator \"op\" should be Hermitian Operator" BosonBath(
        Qobj([0 1; 0 0]),
        [0, 0],
        [0, 0],
        [0, 0],
        [0, 0],
    )
    ################################################

    ################################################
    # Boson bath (RWA)
    b = BosonBathRWA(op, η1, γ3, η2, γ4)
    types = ["bA", "bA", "bA", "bA", "bA", "bE", "bE", "bE", "bE", "bE"]
    @test length(bs) == 10
    cp, cm = correlation_function(b, 0.183183)
    @test_throws ErrorException C(b, [0.183183])
    @test cp[1] ≈ 22.384506765987076 + 0.7399082821797519im
    @test cm[1] ≈ 26.994911851482776 - 0.799138487523946im

    ## check for η and γ list, and coupling operator
    η = Vector{ComplexF64}(undef, length(b))
    γ = Vector{ComplexF64}(undef, length(b))
    for (i, e) in enumerate(b)
        η[i] = e.η
        γ[i] = e.γ
        @test e.op == op
        @test e.types == types[i]
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
    @test show(devnull, MIME("text/plain"), b) === nothing

    ## check for exponents
    @test show(devnull, MIME("text/plain"), b[1]) === nothing
    @test show(devnull, MIME("text/plain"), b[:]) === nothing
    @test show(devnull, MIME("text/plain"), b[1:end]) === nothing

    ## check exceptions
    @test_throws ErrorException BosonBathRWA(op, [0], [0, 0], [0, 0], [0, 0])
    @test_throws ErrorException BosonBathRWA(op, [0, 0], [0], [0, 0], [0, 0])
    @test_throws ErrorException BosonBathRWA(op, [0, 0], [0, 0], [0], [0, 0])
    @test_throws ErrorException BosonBathRWA(op, [0, 0], [0, 0], [0, 0], [0])
    @test_throws ErrorException BosonBathRWA(Qobj([0, 0]), [0, 0], [0, 0], [0, 0], [0, 0])
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
    types = ["fA", "fA", "fA", "fA", "fA", "fE", "fE", "fE", "fE", "fE"]
    @test length(b) == 10
    cp, cm = correlation_function(b, 0.183183)
    @test_throws ErrorException C(b, [0.183183])
    @test cp[1] ≈ 22.384506765987076 + 0.7399082821797519im
    @test cm[1] ≈ 26.994911851482776 - 0.799138487523946im
    for (i, e) in enumerate(b)
        @test e.op == op
        @test e.types == types[i]
    end
    @test show(devnull, MIME("text/plain"), b) === nothing

    ## check for exponents
    @test show(devnull, MIME("text/plain"), b[1]) === nothing
    @test show(devnull, MIME("text/plain"), b[:]) === nothing
    @test show(devnull, MIME("text/plain"), b[1:end]) === nothing

    ## check exceptions
    @test_throws ErrorException FermionBath(op, [0], [0, 0], [0, 0], [0, 0])
    @test_throws ErrorException FermionBath(op, [0, 0], [0], [0, 0], [0, 0])
    @test_throws ErrorException FermionBath(op, [0, 0], [0, 0], [0], [0, 0])
    @test_throws ErrorException FermionBath(op, [0, 0], [0, 0], [0, 0], [0])
    @test_throws ErrorException FermionBath(Qobj([0, 0]), [0, 0], [0, 0], [0, 0], [0, 0])
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
    @test correlation_function(BosonBath(op, η0, γ0), 0.123)[1] ≈ η0[1] * exp(-γ0[1] * 0.123)

    ## check exceptions
    @test_throws ErrorException b[11]
    @test_throws ErrorException b[1:11]
    @test_throws ErrorException b[0:10]
    ################################################
end

# prepare coupling operator and coefficients of exponential-exponential-expansion terms
η1 = [1, 3, 5, 7, 9]
η2 = [2, 4, 6, 8, 10]
γ1 = [0.1, 0.3, 0.5, 0.3, 0.7]
γ2 = [0.1, 0.2, 0.5, 0.6, 0.9]
op = [0 0; 0 0]

################################################
# Boson bath
b = BosonBath(op, η1, η2, combine=false)
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
@test γ == [0.1 + 0.0im,  0.3 + 0.0im, 0.5 + 0.0im, 0.7 + 0.0im]

## check for real and image sperarate case
b = BosonBath(op, η1, γ1, η2, γ2, combine = false)
@test length(b) == 10

## check for combine
b = BosonBath(op, η1, γ1, η2, γ2)
@test length(b) == 7

## check for η and γ list, and coupling operator
η = []
γ = []
for e in b
    push!(η, e.η)
    push!(γ, e.γ)
    @test e.op == [0 0; 0 0]
end
@test η == [10.0 + 0.0im, 9.0 + 0.0im, 4.0 + 0.0im, 8.0 + 0.0im, 10.0 + 0.0im, 1.0 + 2.0im, 5.0 + 6.0im]
@test γ == [ 0.3 + 0.0im, 0.7 + 0.0im, 0.2 + 0.0im, 0.6 + 0.0im,  0.9 + 0.0im, 0.1 + 0.0im, 0.5 + 0.0im]
@test show(devnull, MIME("text/plain"), b) == nothing

## check for exponents
@test show(devnull, MIME("text/plain"), b[1]) == nothing
@test show(devnull, MIME("text/plain"), b[:]) == nothing
@test show(devnull, MIME("text/plain"), b[1:end]) == nothing

## check exceptions
@test_throws ErrorException BosonBath(op, [0], [0, 0])
@test_throws ErrorException BosonBath(op, [0, 0], [0, 0],    [0], [0, 0])
@test_throws ErrorException BosonBath(op, [0, 0],    [0], [0, 0], [0, 0])
@test_throws ErrorException @test_warn "Heom doesn't support matrix type : Vector{Int64}" BosonBath([0, 0], [0, 0], [0, 0], [0, 0], [0, 0])
################################################

################################################
# Fermion bath
b = FermionBath(op, η1, γ1, η2, γ2)
@test length(b) == 10
for e in b
    @test e.op == [0 0; 0 0]
end
@test show(devnull, MIME("text/plain"), b) == nothing

## check for exponents
@test show(devnull, MIME("text/plain"), b[1]) == nothing
@test show(devnull, MIME("text/plain"), b[:]) == nothing
@test show(devnull, MIME("text/plain"), b[1:end]) == nothing

## check exceptions
@test_throws ErrorException FermionBath(op,    [0], [0, 0], [0, 0], [0, 0])
@test_throws ErrorException FermionBath(op, [0, 0],    [0], [0, 0], [0, 0])
@test_throws ErrorException FermionBath(op, [0, 0], [0, 0],    [0], [0, 0])
@test_throws ErrorException FermionBath(op, [0, 0], [0, 0], [0, 0],    [0])
@test_throws ErrorException @test_warn "The size of input matrix should be squared matrix." FermionBath([0, 0], [0, 0], [0, 0], [0, 0], [0, 0])
################################################

################################################
# Exponent
## check exceptions
@test_throws ErrorException b[11]
@test_throws ErrorException b[1:11]
@test_throws ErrorException b[0:10]
################################################
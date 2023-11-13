# Test Schrodinger type HEOM Liouvillian superoperator matrix
t = 10
Hsys = [0 1; 1 0]
L = M_S(Hsys; verbose=false)
@test show(devnull, MIME("text/plain"), L) == nothing
@test size(L) == (4, 4)
@test L.N  == 1
@test nnz(L.data) == 8
ados_list = evolution(L, [1 0; 0 0], 0:1:t; reltol=1e-8, abstol=1e-10, verbose=false)
ados = ados_list[end]
@test ados.dim == L.dim
@test length(ados) == L.N
@test eltype(L) == eltype(ados)
ρ0 = ados[1]
@test getRho(ados) == ρ0
ρ1 = [
          cos(t)^2 -0.5im*sin(2*t);
    0.5im*sin(2*t)        sin(t)^2
]
@test _is_Matrix_approx(ρ0, ρ1)

L =   addBosonDissipator(L, √(0.01) * [1 0; 0 -1])
L = addFermionDissipator(L, √(0.01) * [1 0; 0 -1])
@test nnz(L.data) == 10
ados_list = evolution(L, [1 0; 0 0], 0:0.5:t; reltol=1e-8, abstol=1e-10, verbose=false)
ados = ados_list[end]
ρ0 = ados[1]
@test getRho(ados) == ρ0
ρ1 = [
    0.6711641119639493+0.0im 0.3735796014268062im;
    -0.3735796014268062im    0.32883588803605024
]
@test _is_Matrix_approx(ρ0, ρ1)

## check exceptions
@test_throws BoundsError L[1, 5]
@test_throws BoundsError L[1:5, 2]
@test_throws ErrorException ados[L.N + 1]
@test_throws ErrorException M_S([0, 0]; verbose=false)
# System Hamiltonian and initial state
Hsys = 0.25 * [1 0; 0 -1] + 0.5 * [0 1; 1 0]
ρ0   = [1 0; 0 0]

# Bath properties:
λ  = 0.1
W  = 0.5
kT = 0.5
N  = 2
Q = [1 0; 0 -1]  # System-bath coupling operator
bath = Boson_DrudeLorentz_Pade(Q, λ, W, kT, N)

# HEOM Liouvillian superoperator matrix
tier = 5
L = M_Boson(Hsys, tier, bath; verbose=false)

ados = SteadyState(L, ρ0; verbose=false, reltol=1e-2, abstol=1e-4)
ρs   = getRho(ados)
O = [1 0.5; 0.5 1]
@test Expect(O, ados) ≈ real(tr(O * ρs))
@test Expect(O, ados, take_real=false) ≈ tr(O * ρs)   

ρ1 = getRho(SteadyState(L; verbose=false))
@test _is_Matrix_approx(ρ1, ρs)

mat  = spzeros(ComplexF64, 2, 2)
mat2 = spzeros(ComplexF64, 3, 3)
bathf = Fermion_Lorentz_Pade(mat, 1, 1, 1, 1, 2)
@test_throws ErrorException SteadyState(M_Fermion(mat, 2, bathf, ODD; verbose=false))
@test_throws ErrorException SteadyState(M_Fermion(mat, 2, bathf, ODD; verbose=false), mat)
@test_throws ErrorException SteadyState(L, mat2)
@test_throws ErrorException SteadyState(L, ADOs(zeros(8), 2))
@test_throws ErrorException SteadyState(L, ADOs(ados.data, ados.N, ODD))
@test_throws ErrorException SteadyState(L, HEOMSuperOp(Q, ODD, L) * ados)
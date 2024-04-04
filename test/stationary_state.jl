@time @testset "Stationary state" begin

# System Hamiltonian and initial state
d    = [0 1; 0 0]
Hsys = d' * d
ρ0   = [1 0; 0 0]

# Bath properties:
Γ  =  1.
W  =  1
kT =  0.025851991
μL =  1.
μR = -1.
N  =  3
baths = [
    Fermion_Lorentz_Pade(d, Γ, μL, W, kT, N), 
    Fermion_Lorentz_Pade(d, Γ, μR, W, kT, N)
]

# HEOM Liouvillian superoperator matrix
tier = 5
L = M_Fermion(Hsys, tier, baths; verbose=false)

ados = SteadyState(L, ρ0; verbose=false)
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
@test_throws ErrorException SteadyState(L, HEOMSuperOp(d, ODD, L) * ados)
end
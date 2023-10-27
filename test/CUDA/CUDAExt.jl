using HierarchicalEOM
using CUDA
using CUDA.CUSPARSE
using LinearSolve

# re-define the bath (make the matrix smaller)
λ  = 0.01
W  = 0.5
kT = 0.5
μ  = 0
N  = 3
tier = 3

# System Hamiltonian
Hsys = [
    0 0;
    0 0
]

# system-bath coupling operator
Qb = [
    0 1; 
    1 0
]
Qf = [
    0 1; 
    0 0
]

# initial state
ρ0 = [
    0 0; 
    0 1
]

Bbath = Boson_DrudeLorentz_Pade(Qb, λ, W, kT, N)
Fbath = Fermion_Lorentz_Pade(Qf, λ, μ, W, kT, N)

# Solving time Evolution
## Schrodinger HEOMLS
L_cpu = M_S(Hsys; verbose=false)
L_gpu = cu(L_cpu)
ados_cpu = evolution(L_cpu, ρ0, [0f0, 10f0]; verbose=false)
ados_gpu = evolution(L_gpu, ρ0, [0f0, 10f0]; verbose=false)
@test _is_Matrix_approx(getRho(ados_cpu[end]), getRho(ados_cpu[end]))

## Boson HEOMLS
L_cpu = M_Boson(Hsys, tier, Bbath; verbose=false)
L_gpu = cu(L_cpu)
ados_cpu = evolution(L_cpu, ρ0, [0f0, 10f0]; verbose=false)
ados_gpu = evolution(L_gpu, ρ0, [0f0, 10f0]; verbose=false)
@test _is_Matrix_approx(getRho(ados_cpu[end]), getRho(ados_cpu[end]))

## Boson Fermion HEOMLS
L_cpu = M_Fermion(Hsys, tier, Fbath; verbose=false)
L_gpu = cu(L_cpu)
ados_cpu = evolution(L_cpu, ρ0, [0f0, 10f0]; verbose=false)
ados_gpu = evolution(L_gpu, ρ0, [0f0, 10f0]; verbose=false)
@test _is_Matrix_approx(getRho(ados_cpu[end]), getRho(ados_cpu[end]))

## Boson Fermion HEOMLS
L_cpu = M_Boson_Fermion(Hsys, tier, tier, Bbath, Fbath; verbose=false)
L_gpu = cu(L_cpu)
tlist = 0f0:1f0:10f0
ados_cpu = evolution(L_cpu, ρ0, tlist; verbose=false)
ados_gpu = evolution(L_gpu, ρ0, tlist; verbose=false)
for i in 1:length(tlist)
    @test _is_Matrix_approx(getRho(ados_cpu[i]), getRho(ados_cpu[i]))
end

# Steady State
CUDA.@time ados_cpu = SteadyState(L_cpu; verbose=false)
CUDA.@time ados_gpu = SteadyState(L_gpu; solver=Krylov_GMRES(), verbose=false)
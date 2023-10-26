using HierarchicalEOM
using CUDA
using CUDA.CUSPARSE

# re-define the bath (make the matrix smaller)
λ  = 0.1450
W  = 0.6464
kT = 0.7414
μ  = 0.8787
N  = 3
tier = 2

# System Hamiltonian
Hsys = [
    0.6969 0.4364;
    0.4364 0.3215
]

# system-bath coupling operator
Q = [
               0.1234 0.1357 + 0.2468im; 
    0.1357 - 0.2468im            0.5678
]

# initial state
ρ0 = [
    1 0; 
    0 0
]

Bbath = Boson_DrudeLorentz_Pade(Q, λ, W, kT, N)
Fbath = Fermion_Lorentz_Pade(Q, λ, μ, W, kT, N)

L_cpu = M_Boson_Fermion(Hsys, tier, tier, Bbath, Fbath; verbose=false)
L_gpu = cu(L_cpu)

println(L_cpu)
CUDA.@time ados_cpu = evolution(L_cpu, ρ0, [0f0, 10f0]; verbose=false)
CUDA.@time ados_gpu = evolution(L_gpu, ρ0, [0f0, 10f0]; verbose=false)
@test _is_Matrix_approx(getRho(ados_cpu[end]), getRho(ados_cpu[end]))
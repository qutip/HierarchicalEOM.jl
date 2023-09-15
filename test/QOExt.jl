# System Hamiltonian and initial state
basis = SpinBasis(1//2)
σx = sigmax(basis)
σz = sigmaz(basis)
σm = sigmam(basis)
I2 = identityoperator(basis)

Hsys = 0.25 * σz + 0.5 * σx
ρ0 = dm(Ket(basis, [1, 0]))

λ  = 0.1
W  = 0.5
kT = 0.5
N  = 2
Q = σz  # System-bath coupling operator
bath = Boson_DrudeLorentz_Pade(Q, λ, W, kT, N)

# SteadyState
ados = SteadyState(L, ρ0; verbose=false)

# Evolution
tier = 5
L = M_Boson(Hsys, tier, bath; verbose=false)
O = I2 + 0.5 * σx
Δt    = 10
steps = 1
tlist = 0:Δt:(Δt * steps)
ados_list = evolution(L, ρ0, Δt, steps; verbose=false)
ados_list = evolution(L, ρ0, tlist; verbose=false)
@test Expect(O, ados_list[end]) ≈ real.(tr.(O * QOoperator(ados_list[end][1], Hsys)))

# Power spectral density
a = sigmam(basis)
Hsys = a' * a
λ  = 1e-4
W  = 2e-1
kT = 0.5
N  = 5
bath = Boson_DrudeLorentz_Matsubara((a' + a), λ, W, kT, N)
tier = 3
L = M_Boson(Hsys, tier, bath; verbose=false)
L = addBosonDissipator(L, 1e-3 * a')
L = addTerminator(L, bath)
ados_s = SteadyState(L; verbose=false)
ωlist = [0.9]
@test spectrum(L, ados_s, a, ωlist; verbose=false)[1] ≈ 0.0008880367286438112

# Density of states
e = -5
U = 10
d_up =  σm ⊗ I2
d_dn = -σz ⊗ σm
iden =  I2 ⊗ I2
H0 = e * (d_up' * d_up + d_dn' * d_dn)
H1 = U * (d_up' * d_up * d_dn' * d_dn)
Hsys = H0 + H1
λ   =  1
μ_l =  1
μ_r = -1
W   = 10
kT  =  0.5
N = 5
fuL = Fermion_Lorentz_Pade(d_up, λ, μ_l, W, kT, N)
fdL = Fermion_Lorentz_Pade(d_dn, λ, μ_l, W, kT, N)
fuR = Fermion_Lorentz_Pade(d_up, λ, μ_r, W, kT, N)
fdR = Fermion_Lorentz_Pade(d_dn, λ, μ_r, W, kT, N)
tier = 2
Le = M_Fermion(Hsys, tier, [fuL, fdL, fuR, fdR]; verbose=false)
Lo = M_Fermion(Hsys, tier, [fuL, fdL, fuR, fdR], :odd; verbose=false)
ados_s = SteadyState(Le; verbose=false)
ωlist = [0]
@test spectrum(Lo, ados_s, d_up, ωlist; verbose=false)[1] ≈ 0.17217519700362002
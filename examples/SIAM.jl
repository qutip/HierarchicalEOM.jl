# # [The single-impurity Anderson model](@id exp-SIAM)

# The investigation of the Kondo effect in single-impurity Anderson model is crucial as it serves both as a valuable testing ground for the theories of the Kondo effect and has the potential to lead to a better understanding of this intrinsic many-body phenomena.

using HierarchicalEOM 
import Plots

# ## Hamiltonian
# We consider a single-level electronic system [which can be populated by a spin-up ($\uparrow$) or spin-down ($\downarrow$) electron] coupled to a fermionic reservoir ($\textrm{f}$). The total Hamiltonian is given by $H_{\textrm{T}}=H_\textrm{s}+H_\textrm{f}+H_\textrm{sf}$, where each terms takes the form
# ```math
# \begin{aligned}
# H_{\textrm{s}}  &= \epsilon \left(d^\dagger_\uparrow d_\uparrow + d^\dagger_\downarrow d_\downarrow \right) + U\left(d^\dagger_\uparrow d_\uparrow d^\dagger_\downarrow d_\downarrow\right),\\
# H_{\textrm{f}}  &=\sum_{\sigma=\uparrow,\downarrow}\sum_{k}\epsilon_{\sigma,k}c_{\sigma,k}^{\dagger}c_{\sigma,k},\\
# H_{\textrm{sf}} &=\sum_{\sigma=\uparrow,\downarrow}\sum_{k}g_{k}c_{\sigma,k}^{\dagger}d_{\sigma} + g_{k}^*d_{\sigma}^{\dagger}c_{\sigma,k}.
# \end{aligned}
# ```
# Here, $d_\uparrow$ $(d_\downarrow)$ annihilates a spin-up (spin-down) electron in the system, $\epsilon$ is the energy of the electron, and $U$ is the Coulomb repulsion energy for double occupation. Furthermore, $c_{\sigma,k}$ $(c_{\sigma,k}^{\dagger})$ annihilates (creates) an electron in the state $k$ (with energy $\epsilon_{\sigma,k}$) of the reservoir.
# 
# Now, we can construct the system Hamiltonian
ϵ = -5
U = 10
σm = [0 1; 0  0] ## σ-
σz = [1 0; 0 -1] ## σz
II = [1 0; 0  1] ## identity matrix

## construct the annihilation operator for both spin-up and spin-down
## (utilize Jordan–Wigner transformation)
d_up = kron(     σm, II)
d_dn = kron(-1 * σz, σm)
Hsys = ϵ * (d_up' * d_up + d_dn' * d_dn) + U * (d_up' * d_up * d_dn' * d_dn)

# ## Construct bath objects
# We assume the fermionic reservoir to have a [Lorentzian-shaped spectral density](@ref doc-Fermion-Lorentz), and we utilize the Padé decomposition. Furthermore, the spectral densities depend on the following physical parameters: 
# - the coupling strength $\Gamma$ between system and reservoirs
# - the band-width $W$
# - the temperature $T$
# - the chemical potential $\mu$
# - the total number of exponentials for the reservoir $2(N + 1)$
Γ = 2
μ = 0
W = 10
T = 0.5
N = 5
bath_up = Fermion_Lorentz_Pade(d_up, Γ, μ, W, T, N)
bath_dn = Fermion_Lorentz_Pade(d_dn, Γ, μ, W, T, N)
bath_list = [bath_up, bath_dn]

# ## Construct HEOMLS matrix
# (see also [HEOMLS Matrix for Fermionic Baths](@ref doc-M_Fermion))
tier = 3
M_even = M_Fermion(Hsys, tier, bath_list)
M_odd  = M_Fermion(Hsys, tier, bath_list, :odd)

# ## Solve stationary state of ADOs
# (see also [Stationary State](@ref doc-Stationary-State))
ados_s = SteadyState(M_even)

# ## Calculate density of states (DOS)
# (see also [Spectrum](@ref doc-Spectrum))
ωlist = -10:1:10
dos = spectrum(M_odd, ados_s, d_up, ωlist)

Plots.plot(ωlist, dos)
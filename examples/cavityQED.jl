# # Cavity QED system

# Cavity quantum electrodynamics (cavity QED) is an important topic for studying the interaction between atoms (or other particles) and light confined in a reflective cavity, under conditions where the quantum nature of photons is significant.

using HierarchicalEOM
using LaTeXStrings
import QuantumOptics, Plots

# ## Hamiltonian
# The Jaynes-Cummings model is a standard model in the realm of cavity QED. It illustrates the interaction between a two-level atom ($\textrm{A}$) and a quantized single-mode within a cavity ($\textrm{c}$).

# ```math
# \begin{aligned}
# H_{\textrm{s}}&=H_{\textrm{A}}+H_{\textrm{c}}+H_{\textrm{int}},\\
# H_{\textrm{A}}&=\frac{\omega_A}{2}\sigma_z,\\
# H_{\textrm{c}}&=\omega_{\textrm{c}} a^\dagger a,\\
# H_{\textrm{int}}&=g (a^\dagger\sigma^-+a\sigma^+),
# \end{aligned}
# ```
# where $\sigma^-$ ($\sigma^+$) is the annihilation (creation) operator of the atom, and $a$ ($a^\dagger$) is the annihilation (creation) operator of the cavity.
#   
# Furthermore, we consider the system is coupled to a bosonic reservoir ($\textrm{b}$). The total Hamiltonian is given by $H_{\textrm{Total}}=H_\textrm{s}+H_\textrm{b}+H_\textrm{sb}$, where $H_\textrm{b}$ and $H_\textrm{sb}$ takes the form
# ```math
# \begin{aligned}
# H_{\textrm{b}}    &=\sum_{k}\omega_{k}b_{k}^{\dagger}b_{k},\\
# H_{\textrm{sb}}   &=(a+a^\dagger)\sum_{k}g_{k}(b_k + b_k^{\dagger}).
# \end{aligned}
# ```
# Here, $H_{\textrm{b}}$ describes a bosonic reservoir where $b_{k}$ $(b_{k}^{\dagger})$ is the bosonic annihilation (creation) operator associated to the $k$th mode (with frequency $\omega_{k}$). Also, $H_{\textrm{sb}}$ illustrates the interaction between the cavity and the bosonic reservoir.

# Now, we can build the system Hamiltonian with the package [`QuantumOptics.jl`](https://qojulia.org) (optional) to construct the operators.

# !!! note "Note"
#     Note that `HierarchicalEOM.jl` only accept standard julia bulit-in types. If you use `QuantumOptics` to construct the operators, remember that the matrix (or vector) in standard type are stored in `.data` field of the objects, i.e., `op.data` where `op` is a `QuantumOptics`-type object. Thus, you should take it as the inputs of `HierarchicalEOM.jl` for objects in `QuantumOptics` package.

# !!! compat "Extension for QuantumOptics.jl"
#     `HierarchicalEOM.jl` provides an extension to support `QuantumOptics`-type object, but this feature requires `Julia 1.9+` and `HierarchicalEOM 0.3+`. See [here](@ref doc-ext-QuantumOptics) for more details.

N_photon = 2 ## We truncate the photon number of the system cavity to 2
ωA = 2
ωc = 2
g  = 0.1

## basis
a_basis = QuantumOptics.FockBasis(N_photon)
b_spin  = QuantumOptics.SpinBasis(1//2)

## operators
a_c  = QuantumOptics.destroy(a_basis)
I_c  = QuantumOptics.identityoperator(a_basis)
σz_A = QuantumOptics.sigmaz(b_spin)
σm_A = QuantumOptics.sigmam(b_spin)
I_A  = QuantumOptics.identityoperator(b_spin)

## operators in tensor-space
a  = QuantumOptics.tensor(a_c, I_A)
σz = QuantumOptics.tensor(I_c, σz_A)
σm = QuantumOptics.tensor(I_c, σm_A)

## Hamiltonian
H_A   = 0.5 * ωA * σz
H_c   =       ωc * a' * a
H_int = g * ( a' * σm + a * σm' )

H_s = H_A + H_c + H_int

## initial state
ket0 = QuantumOptics.tensor(QuantumOptics.Ket(a_basis,[1,0,0]),QuantumOptics.Ket(b_spin,[1,0]))
ρ0   = QuantumOptics.dm(ket0);

# ## Construct bath objects
# We assume the bosonic reservoir to have a [Drude-Lorentz Spectral Density](@ref Boson-Drude-Lorentz), and we utilize the Padé decomposition. Furthermore, the spectral densities depend on the following physical parameters: 
# - the coupling strength $\Gamma$ between system and reservoir
# - the band-width $W$
# - the product of the Boltzmann constant $k$ and the absolute temperature $T$ : $kT$
# - the total number of exponentials for the reservoir $(N + 1)$
Γ  = 0.01
W  = 1
kT = 0.025
N  = 20
Bath = Boson_DrudeLorentz_Pade((a + a').data, Γ, W, kT, N)

# Before incorporating the correlation function into the HEOMLS matrix, it is essential to verify if the total number of exponentials for the reservoir sufficiently describes the practical situation.

tlist_test = 0:0.1:10;

Bath_test = Boson_DrudeLorentz_Pade((a + a').data, Γ, W, kT, 1000);
Ct  = C(Bath, tlist_test);
Ct2 = C(Bath_test, tlist_test);

Plots.plot(tlist_test,  real(Ct),  label="N=20 (real part )",  linestyle=:dash,  linewidth=3)
Plots.plot!(tlist_test, real(Ct2), label="N=1000 (real part)", linestyle=:solid, linewidth=3)
Plots.plot!(tlist_test, imag(Ct),  label="N=20 (imag part)",   linestyle=:dash,  linewidth=3)
Plots.plot!(tlist_test, imag(Ct2), label="N=1000 (imag part)", linestyle=:solid, linewidth=3)

Plots.xaxis!("t")
Plots.yaxis!("C(t)")

# ## Construct HEOMLS matrix
# (see also [HEOMLS Matrix for Bosonic Baths](@ref doc-M_Boson))
# Here, we consider an incoherent pumping to the atom, which can be described by an Lindblad dissipator (see [here](@ref doc-Master-Equation) for more details).  
#   
# Furthermore, we set the [important threshold](@ref doc-Importance-Value-and-Threshold) to be `1e-6`.
pump   = 0.01
J_pump = sqrt(pump) * (σm').data

tier = 2
M_Heom = M_Boson(H_s.data, tier, threshold=1e-6, Bath)
M_Heom = addBosonDissipator(M_Heom, J_pump)

# ## Solve time evolution of ADOs
# (see also [Time Evolution](@ref doc-Time-Evolution))
t_list = 0:1:500
evo_H  = evolution(M_Heom, ρ0.data, t_list);

# ## Solve stationary state of ADOs
# (see also [Stationary State](@ref doc-Stationary-State))
steady_H = SteadyState(M_Heom);

# ## Expectation values
# observable of atom: $\sigma_z$
σz_evo_H    = Expect(σz.data, evo_H)
σz_steady_H = Expect(σz.data, steady_H)

# observable of cavity: $a^\dagger a$ (average photon number)
np_evo_H    = Expect((a' * a).data, evo_H)
np_steady_H = Expect((a' * a).data, steady_H)

p1 = Plots.plot(
    t_list,
    [
        σz_evo_H, 
        ones(length(t_list)).* σz_steady_H
    ],
    label=[L"\langle \sigma_z \rangle"  L"\langle \sigma_z \rangle ~~(\textrm{steady})"],
    linewidth=3,
    linestyle=[:solid  :dash]
)
p2 = Plots.plot(
    t_list,
    [
        np_evo_H, 
        ones(length(t_list)).* np_steady_H
    ],
    label=[L"\langle a^\dagger a \rangle"  L"\langle a^\dagger a \rangle ~~(\textrm{steady})"],
    linewidth=3,
    linestyle=[:solid  :dash]
)
Plots.plot(p1, p2, layout=[1,1])
Plots.xaxis!("t")

# ## Power spectral density
# (see also [Spectrum](@ref doc-Spectrum))
ω_list = 1:0.01:3
psd_H = spectrum(M_Heom, steady_H, a.data, ω_list)

Plots.plot(ω_list, psd_H, linewidth=3)
Plots.xaxis!(L"\omega")

# ## Compare with Master Eq. approach
# (see also [HEOMLS for Master Equations](@ref doc-Master-Equation))
#   
# The Lindblad master equations which describs the cavity couples to an extra bosonic reservoir with [Drude-Lorentzian spectrum density](@ref Boson-Drude-Lorentz) is given by

## Drude_Lorentzian spectral density
Drude_Lorentz(ω, Γ, W) = 4 * Γ * W * ω / ( (ω)^2 + (W)^2 )

## Bose-Einstein distribution
n_b(ω, kT) = 1 / (exp(ω / kT) - 1)

## build the jump operators
jump_op = [
    sqrt(Drude_Lorentz(ωc, Γ, W) * (n_b(ωc, kT) + 1)) * a.data,
    sqrt(Drude_Lorentz(ωc, Γ, W) * (n_b(ωc, kT)))     * (a').data,
    J_pump
];

## construct the HEOMLS matrix for master equation
M_master = M_S(H_s.data)
M_master = addBosonDissipator(M_master, jump_op)

## time evolution
evo_M = evolution(M_master, ρ0.data, t_list);

## steady
steady_M = SteadyState(M_master);

## expectation value of σz
σz_evo_M    = Expect(σz.data, evo_M)
σz_steady_M = Expect(σz.data, steady_M)

## average photon number
np_evo_M    = Expect((a' * a).data, evo_M)
np_steady_M = Expect((a' * a).data, steady_M)

p1 = Plots.plot(
    t_list,
    [
        σz_evo_M, 
        ones(length(t_list)).* σz_steady_M
    ],
    label=[L"\langle \sigma_z \rangle"  L"\langle \sigma_z \rangle ~~(\textrm{steady})"],
    linewidth=3,
    linestyle=[:solid  :dash],
)
p2 = Plots.plot(
    t_list,
    [
        np_evo_M, 
        ones(length(t_list)).* np_steady_M
    ],
    label=[L"\langle a^\dagger a \rangle"  L"\langle a^\dagger a \rangle ~~(\textrm{steady})"],
    linewidth=3,
    linestyle=[:solid  :dash],
)
Plots.plot(p1, p2, layout=[1,1])
Plots.xaxis!("t")

# We can also calculate the power spectral density

ω_list = 1:0.01:3
psd_M = spectrum(M_master, steady_M, a.data, ω_list)

Plots.plot(ω_list, psd_M, linewidth=3)
Plots.xaxis!(L"\omega")

# Due to the weak coupling between the system and an extra bosonic environment, the Master equation's outcome is expected to be similar to the results obtained from the HEOM method.
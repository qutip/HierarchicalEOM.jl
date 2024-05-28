# # [Driven Systems and Dynamical Decoupling](@id exp-dynamical-decoupling)
# In this page, we show how to solve the time evolution with time-dependent Hamiltonian problems in hierarchical equations of motion approach.

using HierarchicalEOM
using LaTeXStrings
import Plots

# Here, we study dynamical decoupling which is a common tool used to undo the dephasing effect from the environment even for finite pulse duration.  
#
# ## Hamiltonian
# We consider a two-level system coupled to a bosonic reservoir ($\textrm{b}$). The total Hamiltonian is given by $H_{\textrm{T}}=H_\textrm{s}+H_\textrm{b}+H_\textrm{sb}$, where each terms takes the form
# ```math
# \begin{aligned}
# H_{\textrm{s}}(t) &= H_0 + H_{\textrm{D}}(t),\\
# H_0               &= \frac{\omega_0}{2} \sigma_z,\\
# H_{\textrm{b}}    &=\sum_{k}\omega_{k}b_{k}^{\dagger}b_{k},\\
# H_{\textrm{sb}}   &=\sigma_z\sum_{k}g_{\alpha,k}(b_k + b_k^{\dagger}).
# \end{aligned}
# ```
# Here, $H_{\textrm{b}}$ describes a bosonic reservoir where $b_{k}$ $(b_{k}^{\dagger})$ is the bosonic annihilation (creation) operator associated to the $k$th mode (with frequency $\omega_{k}$).
#   
# Furthermore, to observe the time evolution of the coherence, we consider the initial state to be
# ```math
# ψ(t=0)=\frac{1}{\sqrt{2}}\left(|0\rangle+|1\rangle\right)
# ```
ω0 = 0.0
σz = [1 0; 0 -1]
σx = [0 1; 1 0]
H0 = 0.5 * ω0 * σz

ρ0 = 0.5 * [1 1; 1 1];

# The time-dependent driving term $H_{\textrm{D}}(t)$ has the form
# ```math
# H_{\textrm{D}}(t) = \sum_{n=1}^N f_n(t) \sigma_x 
# ```
# where the pulse is chosen to have duration $\tau$ together with a delay $\Delta$ between each pulses, namely
# ```math
# f_n(t)
# = \begin{cases}
#   V & \textrm{if}~~(n-1)\tau + n\Delta \leq t \leq n (\tau + \Delta),\\
#   0 & \textrm{otherwise}.
#   \end{cases} 
# ```
# Here, we set the period of the pulses to be $\tau V = \pi/2$. Therefore, we consider two scenarios with fast and slow pulses:

## a function which returns the amplitude of the pulse at time t
function pulse(V, Δ, t)
    τ = 0.5 * π / V
    period = τ + Δ

    if (t % period) < τ
        return V
    else
        return 0
    end
end

tlist = 0:0.4:400
amp_fast = 0.50
amp_slow = 0.01
delay = 20

Plots.plot(
    tlist,
    [[pulse(amp_fast, delay, t) for t in tlist], [pulse(amp_slow, delay, t) for t in tlist]],
    label = ["Fast Pulse" "Slow Pulse"],
    linestyle = [:solid :dash],
)

# ## Construct bath objects
# We assume the bosonic reservoir to have a [Drude-Lorentz Spectral Density](@ref Boson-Drude-Lorentz), and we utilize the Padé decomposition. Furthermore, the spectral densities depend on the following physical parameters: 
# - the coupling strength $\Gamma$ between system and reservoir
# - the band-width $W$
# - the product of the Boltzmann constant $k$ and the absolute temperature $T$ : $kT$
# - the total number of exponentials for the reservoir $(N + 1)$
Γ = 0.0005
W = 0.005
kT = 0.05
N = 3
bath = Boson_DrudeLorentz_Pade(σz, Γ, W, kT, N)

# ## Construct HEOMLS matrix
# (see also [HEOMLS Matrix for Bosonic Baths](@ref doc-M_Boson))
# !!! note "Note"
#     Only provide the **time-independent** part of system Hamiltonian when constructing HEOMLS matrices (the time-dependent part should be given when solving the time evolution).
tier = 6
M = M_Boson(H0, tier, bath)

# ## time evolution with time-independent Hamiltonian
# (see also [Time Evolution](@ref doc-Time-Evolution))
noPulseResult = evolution(M, ρ0, tlist);

# ## Solve time evolution with time-dependent Hamiltonian
# (see also [Time Evolution](@ref doc-Time-Evolution))
#   
# We need to provide a user-defined function (named as `H_D` in this case), which must be in the form `H_D(p::Tuple, t)` and returns the time-dependent part of system Hamiltonian (in `AbstractMatrix` type) at any given time point `t`. The parameter `p` should be a `Tuple` which contains all the extra parameters [`V` (amplitude), `Δ` (delay), and `σx` (operator) in this case] for the function `H_D`:
function H_D(p::Tuple, t)
    V, Δ, σx = p
    return pulse(V, Δ, t) * σx
end;

# The parameter tuple `p` will be passed to your function `H_D` directly from the **last required** parameter in `evolution`:
fastTuple = (amp_fast, delay, σx)
slowTuple = (amp_slow, delay, σx)

fastPulseResult = evolution(M, ρ0, tlist, H_D, fastTuple);
slowPulseResult = evolution(M, ρ0, tlist, H_D, slowTuple);

# ## Measure the coherence
# One can use the built-in function [`Expect`](@ref) to calculate the expectation value from a given observable and `ADOs`:

## Define the operator that measures the 0, 1 element of density matrix
ρ01 = [0 1; 0 0]

Plots.plot(
    tlist,
    [Expect(ρ01, fastPulseResult), Expect(ρ01, slowPulseResult), Expect(ρ01, noPulseResult)],
    label = ["Fast Pulse" "Slow Pulse" "no Pulse"],
    linestyle = [:solid :dot :dash],
    linewidth = 3,
    xlabel = L"t",
    ylabel = L"\rho_{01}",
    grid = false,
)

# This example is from QuTiP-BoFiN paper : [Phys. Rev. Research 5, 013181 (2023)](https://link.aps.org/doi/10.1103/PhysRevResearch.5.013181).

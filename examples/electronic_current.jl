# # Electronic Current
# 
# In this example, we demonstrate how to compute an environmental observable: the electronic current. 

using HierarchicalEOM
import LinearAlgebra: tr
import Plots

# ## Hamiltonian
# We consider a single-level charge system coupled to two [left (L) and right (R)] fermionic reservoirs ($\textrm{f}$). The total Hamiltonian is given by $H_{\textrm{T}}=H_\textrm{s}+H_\textrm{f}+H_\textrm{sf}$, where each terms takes the form
# ```math
# H_{\textrm{s}} = \epsilon d^\dagger d,
# ```
# ```math
# H_{\textrm{f}} =\sum_{\alpha=\textrm{L},\textrm{R}}\sum_{k}\epsilon_{\alpha,k}c_{\alpha,k}^{\dagger}c_{\alpha,k},
# ```
# ```math
# H_{\textrm{sf}} =\sum_{\alpha=\textrm{L},\textrm{R}}\sum_{k}g_{\alpha,k}c_{\alpha,k}^{\dagger}d + g_{\alpha,k}^* d^{\dagger}c_{\alpha,k}.
# ```
# Here, $d$ $(d^\dagger)$ annihilates (creates) an electron in the system and $\epsilon$ is the energy of the electron. Furthermore, $c_{\alpha,k}$ $(c_{\alpha,k}^{\dagger})$ annihilates (creates) an electron in the state $k$ (with energy $\epsilon_{\alpha,k}$) of the $\alpha$-th reservoir.
# 
# Now, we can construct the system Hamiltonian

d = [0 1; 0 0] ## annihilation operator of the system electron

## The system Hamiltonian
ϵ = 1. # site energy
Hsys = ϵ * d' * d

## System initial state
ρ0 = [1 0; 0 0];

# ## Construct bath objects
# We assume the fermionic reservoir to have a [Lorentzian-shaped spectral density](@ref doc-Fermion-Lorentz), and we utilize the Padé decomposition. Furthermore, the spectral densities depend on the following physical parameters: 
# - the coupling strength $\Gamma$ between system and reservoirs
# - the band-width $W$
# - the temperature $T$
# - the chemical potential $\mu$
# - the total number of exponentials for the reservoir $2(N + 1)$
Γ = 0.01
W = 1
T = 0.025851991
μL =  1.  # Left  bath
μR = -1.  # Right bath
N = 2
bath_L = Fermion_Lorentz_Pade(d, Γ, μL, W, T, N)
bath_R = Fermion_Lorentz_Pade(d, Γ, μR, W, T, N)
baths = [bath_L, bath_R]

# ## Construct HEOMLS matrix
tier = 5
M = M_Fermion(Hsys, tier, baths)

# ## Solve time evolution of ADOs
tlist = 0:0.5:100
ados_evolution = evolution(M, ρ0, tlist);

# ## Solve stationary state of ADOs
ados_steady = SteadyState(M);

# ## Calculate current
# Within the influence functional approach, the expectation value of the electronic current from the $\alpha$-fermionic bath into the system can be written in terms of the first-level-fermionic ($n=1$) auxiliary density operators, namely
# ```math
# \langle I_\alpha(t) \rangle =(-e) \frac{d\langle \mathcal{N}_\alpha\rangle}{dt}=i e \sum_{q}(-1)^{\delta_{\nu,-}} ~\textrm{Tr}\left[d^{\bar{\nu}}\rho^{(0,1,+)}_{\vert \textbf{q}}(t)\right],
# ```
# where $e$ represents the value of the elementary charge, and $\mathcal{N}_\alpha=\sum_k c^\dagger_{\alpha,k}c_{\alpha,k}$ is the occupation number operator for the $\alpha$-fermionic bath.
# 
# Given an ADOs, we provide a function which calculates the current from the $\alpha$-fermionic bath into the system with the help of [Hierarchy Dictionary](@ref doc-Hierarchy-Dictionary).  
# 
# `bathIdx`:  
# - 1 means 1st bath (bath_L)
# - 2 means 2nd bath (bath_R)
function Ic(ados, M::M_Fermion, bathIdx::Int)
    ## the hierarchy dictionary
    HDict = M.hierarchy

    ## we need all the indices of ADOs for the first level
    idx_list = HDict.lvl2idx[1]
    I = 0.0im
    for idx in idx_list
        ρ1 = ados[idx]  ## 1st-level ADO

        ## find the corresponding bath index (α) and exponent term index (k)
        nvec = HDict.idx2nvec[idx]
        for (α, k, _) in getIndexEnsemble(nvec, HDict.bathPtr)
            if α == bathIdx
                exponent = M.bath[α][k]
                if exponent.types == "fA"     ## fermion-absorption
                    I += tr(exponent.op' * ρ1)
                elseif exponent.types == "fE" ## fermion-emission
                    I -= tr(exponent.op' * ρ1)
                end
                break
            end
        end
    end
    
    e = 1.60218e-19
    ħ = 6.62607015e−34 / (2 * π)
    eV_to_Joule = 1.60218e-19  ## unit conversion
    
    ## (e / ħ) * I  [change unit to μA] 
    return (e / ħ) * real(1im * I) * eV_to_Joule * 1e6
end

## steady current
Is_L = ones(length(tlist)) .* Ic(ados_steady, M, 1)
Is_R = ones(length(tlist)) .* Ic(ados_steady, M, 2)

## time evolution current
Ie_L = []
Ie_R = []
for ados in ados_evolution
    push!(Ie_L, Ic(ados, M, 1))
    push!(Ie_R, Ic(ados, M, 2))
end

Plots.plot(
    tlist, 
    [Ie_L, Ie_R, Is_L, Is_R],
    label=["Bath L" "Bath R" "Bath L (Steady State)" "Bath R (Steady State)"],
    linecolor=[ :blue   :red :blue  :red],
    linestyle=[:solid :solid :dash :dash],
    linewidth=3,
    xlabel="time",
    ylabel="Current",
    grid=false
)

# Note that this example can also be found in [qutip documentation](https://qutip.org/docs/latest/guide/heom/fermionic.html#steady-state-currents)
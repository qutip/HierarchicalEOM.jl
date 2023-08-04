# [Fermionic Bath](@id doc-Fermionic-Bath)

## [Introduction](@id Fermionic-Bath-Introduction)
The [`FermionBath`](@ref) object describes the interaction between the system (``s``) and a exterior fermionic environment (``f``), which can be modeled by
```math
H_{sf}=\sum_k g_k c_k^\dagger d_s + g_k^* c_k d_s^\dagger,
```
where ``g_k`` is the coupling strength and ``c_k (c_k^\dagger)`` annihilates (creates) a fermion in the ``k``-th state of the fermionic environment. Here, ``d_s`` refers to the system-interaction operator and should be an odd-parity operator destroying a fermion in the system.

The effects of a fermionic environment (initially in thermal equilibrium and linearly coupled to the system) are completely encoded in the two-time correlation functions, namely
```math
C^{\nu}(t_{1},t_{2})
=\frac{1}{2\pi}\int_{-\infty}^{\infty} d\omega 
J(\omega)\left[\frac{1-\nu}{2}+\nu n(\omega)
\right]e^{\nu i\omega (t_{1}-t_{2})}.
```
where ``J(\omega)`` is the spectral density of the bath and ``n(\omega)=\{\exp[(\omega-\mu)/k_B T]+1\}^{-1}`` represents the Fermi-Dirac distribution (with chemical potential ``\mu``). Here, ``\nu=+`` and ``\nu=-`` denotes the absorption and emission process of the fermionic system, respectively.

A more practical representation can be found by expressing the correlation function as a sum of exponential terms ([`Exponent`](@ref)), namely
```math
C^{\nu}(t_1, t_2)=\sum_i \eta_i^{\nu} e^{-\gamma_i^{\nu} (t_1-t_2)}.
```
This allows us to define an iterative procedure which leads to the hierarchical equations of motion (HEOM).

## Methods
### Construct FermionBath
One can construct the [`FermionBath`](@ref) object with the coupling operator `ds::AbstractMatrix` and the four lists `η_absorb::AbstractVector`, `γ_absorb::AbstractVector`, `η_emit::AbstractVector` and `γ_emit::AbstractVector` which correspond to the exponential terms ``\{\eta_i^{+}\}_i``, ``\{\gamma_i^{+}\}_i``, ``\{\eta_i^{-}\}_i`` and ``\{\gamma_i^{-}\}_i``, respectively. 
```julia
bath = FermionBath(ds, η_absorb, γ_absorb, η_emit, γ_emit)
```
!!! warning "Warning"
    Here, the length of the four lists (`η_absorb`, `γ_absorb`, `η_emit` and `γ_emit`) should all be the same.

!!! note "Note"
    Instead of analytically solving the correlation function ``C^{\nu=\pm}(t_1, t_2)`` to obtain a sum of exponential terms, one can also use the built-in functions (for different spectral densities ``J(\omega)`` and spectral decomposition methods, which have been analytically solved by the developers already) listed in the end of this page. 

### Print Fermionic Bath
One can check the information of the [`FermionBath`](@ref) by the `print` function, for example:
```julia
print(bath)
```
```
FermionBath object with (system) dim = 2 and 4 exponential-expansion terms
```
Note that [`FermionBath`](@ref) always have even number of exponential terms (half for ``C^{\nu=+}`` and half for ``C^{\nu=-}``)

### Calculate the correlation function
To check whether the exponential terms in the [`FermionBath`](@ref) is correct or not, one can call [`C(bath::FermionBath, tlist::AbstractVector)`](@ref) to calculate the correlation function ``C(t)``, where ``t=t_1-t_2``:
```julia
cp_list, cm_list = C(bath, tlist)
```
Here, `cp_list` and `cm_list` are the lists which contain the value of ``C^{\nu=+}(t)`` and ``C^{\nu=-}(t)`` correspond to the given time series `tlist`, respectively.

### Methods for Exponent
`HierarchicalEOM.jl` also supports users to access the specific exponential term with brakets `[]`. This returns an [`Exponent`](@ref) object, which contains the corresponding value of ``\eta_i^\nu`` and ``\gamma_i^\nu``:
```julia
e = bath[2] # the 2nd-term
print(e)
```
```
Bath Exponent with types = "fA", operator size = (2, 2), η = 0.0 + 3.4090909090909113e-6im, γ = 0.1732050807568877 - 0.005im.
```

One can even obtain the [`Exponent`](@ref) with iterative method:
```julia
for e in bath
    println(e)
end
```
```
Bath Exponent with types = "fA", operator size = (2, 2), η = 6.25e-6 - 3.4090909090909113e-6im, γ = 0.05 - 0.005im.

Bath Exponent with types = "fA", operator size = (2, 2), η = 0.0 + 3.4090909090909113e-6im, γ = 0.1732050807568877 - 0.005im.

Bath Exponent with types = "fE", operator size = (2, 2), η = 6.25e-6 - 3.4090909090909113e-6im, γ = 0.05 + 0.005im.

Bath Exponent with types = "fE", operator size = (2, 2), η = 0.0 + 3.4090909090909113e-6im, γ = 0.1732050807568877 + 0.005im.
```

### Types of Exponent
The different types of the (fermionic-bath) [`Exponent`](@ref):
 - `"fA"` : from absorption fermionic correlation function ``C^{\nu=+}(t_1, t_2)``
 - `"fE"` : from emission fermionic correlation function ``C^{\nu=-}(t_1, t_2)``

## [Lorentz Spectral Density](@id doc-Fermion-Lorentz)
```math
J(\omega)=\frac{\Gamma W^2}{(\omega-\mu)^2+W^2}
```
Here, ``\Gamma`` represents the coupling strength between system and the fermionic environment with chemical potential ``\mu`` and band-width ``W``.

### Matsubara Expansion
With Matsubara Expansion, the correlation function can be analytically solved and expressed as follows:
```math
C^{\nu}(t_1, t_2)=\sum_{l=1}^{\infty} \eta_l^\nu \exp(-\gamma_l^\nu (t_1-t_2))
```
with
```math
\begin{aligned}
\gamma_{1}^{\nu} &= W-\nu i \mu,\\
\eta_{1}^{\nu} &= \frac{\Gamma W}{2} f\left(\frac{iW}{k_B T}\right),\\
\gamma_{l\neq 1}^{\nu} &= \zeta_l k_B T - \nu i \mu,\\
\eta_{l\neq 1}^{\nu} &= -i k_B T \cdot \frac{\Gamma W^2}{-(\zeta_l k_B T)^2+W^2},\\
f(x) &= \{\exp(x) + 1\}^{-1},
\end{aligned}
```
where ``\zeta_l=(2 l - 1)\pi``. This can be constructed by the built-in function [`Fermion_Lorentz_Matsubara`](@ref):
```julia
ds # coupling operator
Γ  # coupling strength
μ  # chemical potential of the environment
W  # band-width  of the environment
T  # temperature of the environment
N  # Number of exponential terms for each correlation functions (C^{+} and C^{-})
bath = Fermion_Lorentz_Matsubara(ds, Γ, μ, W, T, N - 1)
```

### Padé Expansion
With Padé Expansion, the correlation function can be analytically solved and expressed as the following exponential terms:
```math
C^\nu(t_1, t_2)=\sum_{l=1}^{\infty} \eta_l^\nu \exp(-\gamma_l^\nu (t_1-t_2))
```
with
```math
\begin{aligned}
\gamma_{1}^{\nu} &= W-\nu i \mu,\\
\eta_{1}^{\nu} &= \frac{\Gamma W}{2} f\left(\frac{iW}{k_B T}\right),\\
\gamma_{l\neq 1}^{\nu} &= \zeta_l k_B T - \nu i \mu,\\
\eta_{l\neq 1}^{\nu} &= -i \kappa_l k_B T \cdot \frac{\Gamma W^2}{-(\zeta_l k_B T)^2+W^2},\\
f(x) &= \frac{1}{2}-\sum_{l=2}^{N} \frac{2\kappa_l x}{x^2+\zeta_l^2},
\end{aligned}
```
where the parameters ``\kappa_l`` and ``\zeta_l`` are described in [J. Chem. Phys. 134, 244106 (2011)](https://doi.org/10.1063/1.3602466) and ``N`` represents the number of exponential terms for ``C^{\nu=\pm}``. This can be constructed by the built-in function [`Fermion_Lorentz_Pade`](@ref):
```julia
ds # coupling operator
Γ  # coupling strength
μ  # chemical potential of the environment
W  # band-width  of the environment
T  # temperature of the environment
N  # Number of exponential terms for each correlation functions (C^{+} and C^{-})
bath = Fermion_Lorentz_Pade(ds, Γ, μ, W, T, N - 1)
```
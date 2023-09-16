# [Lorentz Spectral Density](@id doc-Fermion-Lorentz)
```math
J(\omega)=\frac{\Gamma W^2}{(\omega-\mu)^2+W^2}
```
Here, ``\Gamma`` represents the coupling strength between system and the fermionic environment with chemical potential ``\mu`` and band-width ``W``.

## Matsubara Expansion
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
kT # the product of the Boltzmann constant k and the absolute temperature T
N  # Number of exponential terms for each correlation functions (C^{+} and C^{-})
bath = Fermion_Lorentz_Matsubara(ds, Γ, μ, W, kT, N - 1)
```

## Padé Expansion
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
kT # the product of the Boltzmann constant k and the absolute temperature T
N  # Number of exponential terms for each correlation functions (C^{+} and C^{-})
bath = Fermion_Lorentz_Pade(ds, Γ, μ, W, kT, N - 1)
```
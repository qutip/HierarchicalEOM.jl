# [Drude-Lorentz Spectral Density](@id Boson-Drude-Lorentz)
```math
J(\omega)=\frac{4\Delta W\omega}{\omega^2+W^2}
```
Here, ``\Delta`` represents the coupling strength between system and the bosonic environment with band-width ``W``.

## Matsubara Expansion
With Matsubara Expansion, the correlation function can be analytically solved and expressed as follows:
```math
C(t_1, t_2)=\sum_{l=1}^{\infty} \eta_l \exp(-\gamma_l (t_1-t_2))
```
with
```math
\begin{aligned}
\gamma_{1} &= W,\\
\eta_{1} &= \Delta W\left[-i+\cot\left(\frac{W}{2 k_B T}\right)\right],\\
\gamma_{l\neq 1} &= 2\pi l k_B T,\\
\eta_{l\neq 1} &= -2 k_B T \cdot \frac{2\Delta W \cdot \gamma_l}{-\gamma_l^2 + W^2}.
\end{aligned}
```
This can be constructed by the built-in function [`Boson_DrudeLorentz_Matsubara`](@ref):
```julia
Vs # coupling operator
Δ  # coupling strength
W  # band-width  of the environment
kT # the product of the Boltzmann constant k and the absolute temperature T
N  # Number of exponential terms
bath = Boson_DrudeLorentz_Matsubara(Vs, Δ, W, kT, N - 1)
```

## Padé Expansion
With Padé Expansion, the correlation function can be analytically solved and expressed as the following exponential terms:
```math
C(t_1, t_2)=\sum_{l=1}^{\infty} \eta_l \exp(-\gamma_l (t_1-t_2))
```
with
```math
\begin{aligned}
\gamma_{1} &= W,\\
\eta_{1} &= \Delta W\left[-i+\cot\left(\frac{W}{2 k_B T}\right)\right],\\
\gamma_{l\neq 1} &= \zeta_l k_B T,\\
\eta_{l\neq 1} &= -2 \kappa_l k_B T \cdot \frac{2\Delta W \cdot \zeta_l k_B T}{-(\zeta_l k_B T)^2 + W^2},
\end{aligned}
```
where the parameters ``\kappa_l`` and ``\zeta_l`` are described in [J. Chem. Phys. 134, 244106 (2011)](https://doi.org/10.1063/1.3602466). This can be constructed by the built-in function [`Boson_DrudeLorentz_Pade`](@ref):
```julia
Vs # coupling operator
Δ  # coupling strength
W  # band-width  of the environment
kT # the product of the Boltzmann constant k and the absolute temperature T
N  # Number of exponential terms
bath = Boson_DrudeLorentz_Pade(Vs, Δ, W, kT, N - 1)
```
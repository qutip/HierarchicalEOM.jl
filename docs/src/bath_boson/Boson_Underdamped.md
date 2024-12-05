# [Underdamped Spectral Density](@id Boson-Underdamped)
```math
J(\omega) = 2 \pi \sum_k |g_k|^2 \delta(\omega-\omega_k) = \frac{2 \Delta^2 W \omega}{(\omega^2 - \omega_0^2)^2 + \omega^2 W^2}
```
Here, ``\Delta`` represents the coupling strength between system and the bosonic environment with band-width ``W`` and resonance frequency ``\omega_0``.

## Matsubara Expansion
With Matsubara Expansion, the correlation function can be analytically solved and expressed as follows:
```math
C(t_1, t_2) = C^\mathrm{R}(t_1, t_2) + iC^\mathrm{I}(t_1, t_2) = \sum_{l=1}^{\infty} \eta_l^\mathrm{R} \exp(-\gamma_l^\mathrm{R} (t_1-t_2)) + \sum_{l=1}^{2} \eta_l^\mathrm{I} \exp(-\gamma_l^\mathrm{I} (t_1-t_2))
```
with
```math
\begin{aligned}
\gamma_{1}^\mathrm{R} &= -i\Omega + \frac{W}{2},\\
\eta_{1}^\mathrm{R} &= \frac{\Delta^2}{4\Omega}\coth\left[\frac{1}{2 k_B T}\left(\Omega + i\frac{W}{2}\right)\right],\\
\gamma_{2}^\mathrm{R} &= i\Omega + \frac{W}{2},\\
\eta_{2}^\mathrm{R} &= \frac{\Delta^2}{4\Omega}\coth\left[\frac{1}{2 k_B T}\left(\Omega - i\frac{W}{2}\right)\right],\\
\gamma_{l}^\mathrm{R} &= 2\pi l k_B T ~~\forall~~ l \geq 3,\\
\eta_{l}^\mathrm{R} &= -2 k_B T \cdot \frac{\Delta^2 W \cdot \gamma_l^\mathrm{R}}{\left[\left(\Omega + i\frac{W}{2}\right)^2 + {\gamma_l^\mathrm{R}}^2\right]\left[\left(\Omega - i\frac{W}{2}\right)^2 + {\gamma_l^\mathrm{R}}^2\right]} ~~\forall~~ l \geq 3,\\
\gamma_{1}^\mathrm{I} &= i\Omega + \frac{W}{2},\\
\eta_{1}^\mathrm{I} &= i\frac{\Delta^2}{4\Omega},\\
\gamma_{2}^\mathrm{I} &= -i\Omega + \frac{W}{2},\\
\eta_{2}^\mathrm{I} &= -i\frac{\Delta^2}{4\Omega},
\end{aligned}
```
where ``\Omega = \sqrt{\omega_0^2 - (W/2)^2}``.
This can be constructed by the built-in function [`Boson_Underdamped_Matsubara`](@ref):
```julia
Vs # coupling operator
Δ  # coupling strength
W  # band-width of the environment
ω0 # resonance frequency of the environment
kT # the product of the Boltzmann constant k and the absolute temperature T
N  # Number of exponential terms
bath = Boson_Underdamped_Matsubara(Vs, Δ, W, ω0, kT, N - 2)
```
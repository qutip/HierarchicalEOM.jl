# [Bosonic Bath](@id doc-Bosonic-Bath)
## [Introduction](@id Bosonic-Bath-Introduction)
The [`BosonBath`](@ref) object describes the interaction between the system (``s``) and a exterior bosonic environment (``b``), which can be modeled by
```math
H_{sb}=V_{s}\sum_k g_k (b_k + b_k^\dagger),
```
in terms of the coupling strength ``g_k`` and the annihilation (creation) operator ``b_k (b_k^\dagger)`` associated to the ``k``-th mode of the bosonic environment. Here, ``V_s`` refers to the system-interaction operator. In particular, ``V_s`` must be a Hermitian operator which can act on **both bosonic and fermionic systems degree of freedom**. In the fermionic system case, ``V_s`` must have even parity to be compatible with charge conservation.

The effects of a bosonic environment (initially in thermal equilibrium and linearly coupled to the system) are completely encoded in the two-time correlation functions, namely
```math
C(t_1, t_2)
=\frac{1}{2\pi}\int_{0}^{\infty} d\omega J(\omega)\left[n(\omega)e^{i\omega (t_1-t_2)}+(n(\omega)+1)e^{-i\omega (t_1-t_2)}\right],
```
where ``J(\omega)`` is the spectral density of the bath and ``n(\omega)=\{\exp(\omega/k_B T)-1\}^{-1}`` represents the Bose-Einstein distribution.

A more practical representation can be found by expressing the correlation function as a sum of exponential terms ([`Exponent`](@ref)), namely
```math
C(t_1, t_2)=\sum_i \eta_i e^{-\gamma_i (t_1-t_2)}.
```
This allows us to define an iterative procedure which leads to the hierarchical equations of motion (HEOM).

## Methods

### Construct BosonBath (with real and imaginary parts are combined)
One can construct the [`BosonBath`](@ref) object with the coupling operator `Vs::AbstractMatrix` and the two lists `η::AbstractVector` and `γ::AbstractVector` which corresponds to the exponential terms ``\{\eta_i\}_i`` and ``\{\gamma_i\}_i``, respectively.
```julia
bath = BosonBath(Vs, η, γ)
```
!!! warning "Warning"
    Here, the length of `η` and `γ` should be the same.

### Construct BosonBath (with real and imaginary parts are separated)
When ``\gamma_i \neq \gamma_i^*``, a closed form for the HEOM can be obtained by further decomposing ``C(t_1, t_2)`` into its real (R) and imaginary (I) parts as
```math
C(t_1, t_2)=\sum_{u=\textrm{R},\textrm{I}}(\delta_{u, \textrm{R}} + i\delta_{u, \textrm{I}})C^{u}(t_1, t_2)
```
where ``\delta`` is the Kronecker delta function and ``C^{u}(t_1, t_2)=\sum_i \eta_i^u \exp(-\gamma_i^u (t_1-t_2))``

In this case, the [`BosonBath`](@ref) object can be constructed by the following method:
```julia
bath = BosonBath(Vs, η_real, γ_real, η_imag, γ_imag)
```
!!! warning "Warning"
    Here, the length of `η_real` and `γ_real` should be the same.  
    Also, the length of `η_imag` and `γ_imag` should be the same.
Here, `η_real::AbstractVector`, `γ_real::AbstractVector`, `η_imag::AbstractVector` and `γ_imag::AbstractVector` correspond to the exponential terms ``\{\eta_i^{\textrm{R}}\}_i``, ``\{\gamma_i^{\textrm{R}}\}_i``, ``\{\eta_i^{\textrm{I}}\}_i`` and ``\{\gamma_i^{\textrm{I}}\}_i``, respectively.

!!! note "Note"
    Instead of analytically solving the correlation function ``C(t_1, t_2)`` to obtain a sum of exponential terms, one can also use the built-in functions (for different spectral densities ``J(\omega)`` and spectral decomposition methods, which have been analytically solved by the developers already) listed in the end of this page. 

### Print Bosonic Bath
One can check the information of the [`BosonBath`](@ref) by the `print` function, for example:
```julia
print(bath)
```
```
BosonBath object with (system) dim = 2 and 4 exponential-expansion terms
```

### Calculate the correlation function
To check whether the exponential terms in the [`BosonBath`](@ref) is correct or not, one can call [`C(bath::BosonBath, tlist::AbstractVector)`](@ref) to calculate the correlation function ``C(t)``, where ``t=t_1-t_2``:
```julia
c_list = C(bath, tlist)
```
Here, `c_list` is a list which contains the value of ``C(t)`` corresponds to the given time series `tlist`.

### Methods for Exponent
`HierarchicalEOM.jl` also supports users to access the specific exponential term with brakets `[]`. This returns an [`Exponent`](@ref) object, which contains the corresponding value of ``\eta_i`` and ``\gamma_i``:
```julia
e = bath[2] # the 2nd-term
print(e)
```
```
Bath Exponent with types = "bRI", operator size = (2, 2), η = 1.5922874021206546e-6 + 0.0im, γ = 0.3141645167860635 + 0.0im.
```

One can even obtain the [`Exponent`](@ref) with iterative method:
```julia
for e in bath
    println(e)
end
```
```
Bath Exponent with types = "bRI", operator size = (2, 2), η = 4.995832638723504e-5 - 2.5e-6im, γ = 0.005 + 0.0im.

Bath Exponent with types = "bRI", operator size = (2, 2), η = 1.5922874021206546e-6 + 0.0im, γ = 0.3141645167860635 + 0.0im.

Bath Exponent with types = "bRI", operator size = (2, 2), η = 1.0039844180003819e-6 + 0.0im, γ = 0.6479143347831898 + 0.0im.

Bath Exponent with types = "bRI", operator size = (2, 2), η = 3.1005439801387293e-6 + 0.0im, γ = 1.8059644711829272 + 0.0im.
```

### Types of Exponent
The different types of the (bosonic-bath) [`Exponent`](@ref):
 - `"bR"` : from real part of bosonic correlation function ``C^{u=\textrm{R}}(t_1, t_2)``
 - `"bI"` : from imaginary part of bosonic correlation function ``C^{u=\textrm{I}}(t_1, t_2)``
 - `"bRI"` : from combined (real and imaginary part) bosonic bath correlation function ``C(t_1, t_2)``

## [Drude-Lorentz Spectral Density](@id Boson-Drude-Lorentz)
```math
J(\omega)=\frac{4\Delta W\omega}{\omega^2+W^2}
```
Here, ``\Delta`` represents the coupling strength between system and the bosonic environment with band-width ``W``.

### Matsubara Expansion
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

### Padé Expansion
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
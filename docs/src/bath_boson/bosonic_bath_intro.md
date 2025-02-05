# [Bosonic Bath](@id doc-Bosonic-Bath)
## [Overview](@id Bosonic-Bath-Overview)
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
where ``J(\omega)=2\pi\Sigma_k |g_k|^2 \delta(\omega-\omega_k)`` is the spectral density of the bath and ``n(\omega)=\{\exp(\omega/k_B T)-1\}^{-1}`` represents the Bose-Einstein distribution.

A more practical representation can be found by expressing the correlation function as a sum of exponential terms ([`Exponent`](@ref)), namely
```math
C(t_1, t_2)=\sum_i \eta_i e^{-\gamma_i (t_1-t_2)}.
```
This allows us to define an iterative procedure which leads to the hierarchical equations of motion (HEOM).

## Construct BosonBath (with real and imaginary parts are combined)
One can construct the [`BosonBath`](@ref) object with the coupling operator `Vs::QuantumObject` and the two lists `η::AbstractVector` and `γ::AbstractVector` which corresponds to the exponential terms ``\{\eta_i\}_i`` and ``\{\gamma_i\}_i``, respectively.
```julia
bath = BosonBath(Vs, η, γ)
```
!!! warning "Warning"
    Here, the length of `η` and `γ` should be the same.

## Construct BosonBath (with real and imaginary parts are separated)
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
    Instead of analytically solving the correlation function ``C(t_1, t_2)`` to obtain a sum of exponential terms, one can also use the built-in functions (for different spectral densities ``J(\omega)`` and spectral decomposition methods, which have been analytically solved by the developers already). See the other categories of the Bosonic Bath in the sidebar for more details.

## Print Bosonic Bath
One can check the information of the [`BosonBath`](@ref) by the `print` function, for example:
```julia
print(bath)
```
```
BosonBath object with 4 exponential-expansion terms
```

## Calculate the correlation function
To check whether the exponential terms in the [`BosonBath`](@ref) is correct or not, one can call [`correlation_function`](@ref) to calculate the correlation function ``C(t)``, where ``t=t_1-t_2``:
```julia
c_list = correlation_function(bath, tlist)
```
Here, `c_list` is a list which contains the value of ``C(t)`` corresponds to the given time series `tlist`.

## Exponent
`HierarchicalEOM.jl` also supports users to access the specific exponential term with brackets `[]`. This returns an [`Exponent`](@ref) object, which contains the corresponding value of ``\eta_i`` and ``\gamma_i``:
```julia
e = bath[2] # the 2nd-term
print(e)
```
```
Bath Exponent with types = "bRI", η = 1.5922874021206546e-6 + 0.0im, γ = 0.3141645167860635 + 0.0im.
```

The different types of the (bosonic-bath) [`Exponent`](@ref):
 - `"bR"` : from real part of bosonic correlation function ``C^{u=\textrm{R}}(t_1, t_2)``
 - `"bI"` : from imaginary part of bosonic correlation function ``C^{u=\textrm{I}}(t_1, t_2)``
 - `"bRI"` : from combined (real and imaginary part) bosonic bath correlation function ``C(t_1, t_2)``

One can even obtain the [`Exponent`](@ref) with iterative method:
```julia
for e in bath
    println(e)
end
```
```
Bath Exponent with types = "bRI", η = 4.995832638723504e-5 - 2.5e-6im, γ = 0.005 + 0.0im.

Bath Exponent with types = "bRI", η = 1.5922874021206546e-6 + 0.0im, γ = 0.3141645167860635 + 0.0im.

Bath Exponent with types = "bRI", η = 1.0039844180003819e-6 + 0.0im, γ = 0.6479143347831898 + 0.0im.

Bath Exponent with types = "bRI", η = 3.1005439801387293e-6 + 0.0im, γ = 1.8059644711829272 + 0.0im.
```
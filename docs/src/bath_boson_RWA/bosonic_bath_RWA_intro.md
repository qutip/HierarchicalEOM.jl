# [Bosonic Bath (under rotating wave approximation)](@id doc-Bosonic-Bath-RWA)

## [Overview](@id Bosonic-Bath-RWA-Overview)
This describes the interaction between the system (``s``) and a exterior bosonic environment (``b``) under the rotating wave approximation (RWA), which can be modeled by
```math
H_{sb}=\sum_k g_k b_k^\dagger a_s + g_k^* b_k a_s^\dagger,
```
where ``g_k`` is the coupling strength and ``b_k (b_k^\dagger)`` is the annihilation (creation) operator for the ``k``-th mode of the bosonic environment. Here, ``a_s`` refers to the system annihilation operator.

The effects of a bosonic environment (initially in thermal equilibrium, linearly coupled to the system, and under the rotating wave approximation) are completely encoded in the two-time correlation functions, namely
```math
C^{\nu}(t_{1},t_{2})
=\frac{1}{2\pi}\int_{0}^{\infty} d\omega 
J(\omega)\left[\delta_{\nu,-1}+ n(\omega)
\right]e^{\nu i\omega (t_{1}-t_{2})}.
```
where ``J(\omega)=2\pi\Sigma_k |g_k|^2 \delta(\omega-\omega_k)`` is the spectral density of the bath and ``n(\omega)=\{\exp[\omega/k_B T]-1\}^{-1}`` represents the Bose-Einstein distribution. Here, ``\nu=+`` and ``\nu=-`` denotes the absorption and emission process of the bosonic system, respectively.

A more practical representation can be found by expressing the correlation function as a sum of exponential terms ([`Exponent`](@ref)), namely
```math
C^{\nu}(t_1, t_2)=\sum_i \eta_i^{\nu} e^{-\gamma_i^{\nu} (t_1-t_2)}.
```
This allows us to define an iterative procedure which leads to the hierarchical equations of motion (HEOM).

## Construct BosonBath (RWA)
One can construct the [`BosonBath`](@ref) object under RWA by calling the function [`BosonBathRWA`](@ref) together with the following parameters: system annihilation operator `a_s::QuantumObject` and the four lists `η_absorb::AbstractVector`, `γ_absorb::AbstractVector`, `η_emit::AbstractVector` and `γ_emit::AbstractVector` which correspond to the exponential terms ``\{\eta_i^{+}\}_i``, ``\{\gamma_i^{+}\}_i``, ``\{\eta_i^{-}\}_i`` and ``\{\gamma_i^{-}\}_i``, respectively. 
```julia
bath = BosonBathRWA(a_s, η_absorb, γ_absorb, η_emit, γ_emit)
```
!!! warning "Warning"
    Here, the length of the four lists (`η_absorb`, `γ_absorb`, `η_emit` and `γ_emit`) should all be the same. Also, all the elements in `γ_absorb` should be complex conjugate of the corresponding elements in `γ_emit`.

## Print Bosonic Bath
One can check the information of the [`BosonBath`](@ref) by the `print` function, for example:
```julia
print(bath)
```
```
BosonBath object with 4 exponential-expansion terms
```
Note that [`BosonBath`](@ref) under RWA always have even number of exponential terms (half for ``C^{\nu=+}`` and half for ``C^{\nu=-}``)

## Calculate the correlation function
To check whether the exponential terms in the [`FermionBath`](@ref) is correct or not, one can call [`correlation_function`](@ref) to calculate the correlation function ``C(t)``, where ``t=t_1-t_2``:
```julia
cp_list, cm_list = correlation_function(bath, tlist)
```
Here, `cp_list` and `cm_list` are the lists which contain the value of ``C^{\nu=+}(t)`` and ``C^{\nu=-}(t)`` correspond to the given time series `tlist`, respectively.

## Exponent
`HierarchicalEOM.jl` also supports users to access the specific exponential term with brackets `[]`. This returns an [`Exponent`](@ref) object, which contains the corresponding value of ``\eta_i^\nu`` and ``\gamma_i^\nu``:
```julia
e = bath[2] # the 2nd-term
print(e)
```
```
Bath Exponent with types = "bA", η = 0.0 + 3.4090909090909113e-6im, γ = 0.1732050807568877 - 0.005im.
```

The different types of the (bosonic-bath under RWA) [`Exponent`](@ref):
 - `"bA"` : from absorption bosonic correlation function ``C^{\nu=+}(t_1, t_2)``
 - `"bE"` : from emission bosonic correlation function ``C^{\nu=-}(t_1, t_2)``

One can even obtain the [`Exponent`](@ref) with iterative method:
```julia
for e in bath
    println(e)
end
```
```
Bath Exponent with types = "bA", η = 6.25e-6 - 3.4090909090909113e-6im, γ = 0.05 - 0.005im.

Bath Exponent with types = "bA", η = 0.0 + 3.4090909090909113e-6im, γ = 0.1732050807568877 - 0.005im.

Bath Exponent with types = "bE", η = 6.25e-6 - 3.4090909090909113e-6im, γ = 0.05 + 0.005im.

Bath Exponent with types = "bE", η = 0.0 + 3.4090909090909113e-6im, γ = 0.1732050807568877 + 0.005im.
```
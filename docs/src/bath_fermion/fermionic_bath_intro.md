# [Fermionic Bath](@id doc-Fermionic-Bath)

## [Overview](@id Fermionic-Bath-Overview)
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
where ``J(\omega)=2\pi\Sigma_k |g_k|^2 \delta(\omega-\omega_k)`` is the spectral density of the bath and ``n(\omega)=\{\exp[(\omega-\mu)/k_B T]+1\}^{-1}`` represents the Fermi-Dirac distribution (with chemical potential ``\mu``). Here, ``\nu=+`` and ``\nu=-`` denotes the absorption and emission process of the fermionic system, respectively.

A more practical representation can be found by expressing the correlation function as a sum of exponential terms ([`Exponent`](@ref)), namely
```math
C^{\nu}(t_1, t_2)=\sum_i \eta_i^{\nu} e^{-\gamma_i^{\nu} (t_1-t_2)}.
```
This allows us to define an iterative procedure which leads to the hierarchical equations of motion (HEOM).

## Construct FermionBath
One can construct the [`FermionBath`](@ref) object with the system annihilation operator `ds::QuantumObject` and the four lists `η_absorb::AbstractVector`, `γ_absorb::AbstractVector`, `η_emit::AbstractVector` and `γ_emit::AbstractVector` which correspond to the exponential terms ``\{\eta_i^{+}\}_i``, ``\{\gamma_i^{+}\}_i``, ``\{\eta_i^{-}\}_i`` and ``\{\gamma_i^{-}\}_i``, respectively. 
```julia
bath = FermionBath(ds, η_absorb, γ_absorb, η_emit, γ_emit)
```
!!! warning "Warning"
    Here, the length of the four lists (`η_absorb`, `γ_absorb`, `η_emit` and `γ_emit`) should all be the same. Also, all the elements in `γ_absorb` should be complex conjugate of the corresponding elements in `γ_emit`.

!!! note "Note"
    Instead of analytically solving the correlation function ``C^{\nu=\pm}(t_1, t_2)`` to obtain a sum of exponential terms, one can also use the built-in functions (for different spectral densities ``J(\omega)`` and spectral decomposition methods, which have been analytically solved by the developers already). See the other categories of the Fermionic Bath in the sidebar for more details.

## Print Fermionic Bath
One can check the information of the [`FermionBath`](@ref) by the `print` function, for example:
```julia
print(bath)
```
```
FermionBath object with 4 exponential-expansion terms
```
Note that [`FermionBath`](@ref) always have even number of exponential terms (half for ``C^{\nu=+}`` and half for ``C^{\nu=-}``)

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
Bath Exponent with types = "fA", η = 0.0 + 3.4090909090909113e-6im, γ = 0.1732050807568877 - 0.005im.
```

The different types of the (fermionic-bath) [`Exponent`](@ref):
 - `"fA"` : from absorption fermionic correlation function ``C^{\nu=+}(t_1, t_2)``
 - `"fE"` : from emission fermionic correlation function ``C^{\nu=-}(t_1, t_2)``

One can even obtain the [`Exponent`](@ref) with iterative method:
```julia
for e in bath
    println(e)
end
```
```
Bath Exponent with types = "fA", η = 6.25e-6 - 3.4090909090909113e-6im, γ = 0.05 - 0.005im.

Bath Exponent with types = "fA", η = 0.0 + 3.4090909090909113e-6im, γ = 0.1732050807568877 - 0.005im.

Bath Exponent with types = "fE", η = 6.25e-6 - 3.4090909090909113e-6im, γ = 0.05 + 0.005im.

Bath Exponent with types = "fE", η = 0.0 + 3.4090909090909113e-6im, γ = 0.1732050807568877 + 0.005im.
```
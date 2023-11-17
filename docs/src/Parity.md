# [Parity Support](@id doc-Parity)
## Introduction
When the system Hamiltonian contains fermionic systems, the HEOMLS matrix ``\hat{\mathcal{M}}`` might be constructed into a different one depend on the parity of the input operator which ``\hat{\mathcal{M}}`` is acting on. This dependence intuitively originates from the properties of partial traces over composite fermionic spaces, where operators do not necessarily commute. 

As an example, for an environment made out of a single fermion, the reduced matrix elements ``\langle{i}|\rho_\textrm{s}^p|{j}\rangle`` (in a basis  labeled by ``\langle i|`` and ``|{j}\rangle``) involve the perturbative sum of expressions of the form ``\langle i| (c \tilde{\rho}_\textrm{e} \tilde{\rho}_\textrm{s}^p c^\dagger+\tilde{\rho}_\textrm{e} \tilde{\rho}_\textrm{s}^p)|{j}\rangle`` (in terms of environmental operators ``\tilde{\rho}_{\textrm{e}}``, system operators ``\tilde{\rho}_\textrm{s}^p`` with parity ``p``, and the environment-annihilation operator ``c``). These quantities depend on the commutator between ``\tilde{\rho}_\textrm{s}^p`` and ``c``, which is trivial only for `EVEN`-parity (``p=+``). In the `ODD`-parity (``p=-``) case, the partial trace over the environment requires further anti-commutations, ultimately resulting in extra minus signs in the expression for the effective propagator describing the reduced dynamics. 

It is important to explicitly note that, here, by `parity` we do not refer to the presence of an odd or even number of fermions in the system but, rather, to the number of fermionic (annihilation or creation) operators needed to represent ``\rho_\textrm{s}^p``. The reduced density matrix of the system should be an `EVEN`-parity operator and can be expressed as ``\rho_{\textrm{s}}^{p=+}(t)``. However, there are some situations (for example, [calculating density of states for fermionic systems](@ref doc-DOS)) where ``\hat{\mathcal{M}}`` is acting on `ODD`-parity ADOs, e.g., ``\rho_{\textrm{s}}^{p=-}(t)=d_{\textrm{s}}\rho_{\textrm{s}}^{+}(t)`` or ``\rho_{\textrm{s}}^{p=-}(t)=d_{\textrm{s}}^\dagger\rho_{\textrm{s}}^{+}(t)``, where ``d_{\textrm{s}}`` is an annihilation operator acting on fermionic systems.

## Parity support for HEOMLS
One can specify the parameter `parity::AbstractParity` in the function of constructing ``\hat{\mathcal{M}}`` which describes the dynamics of [`EVEN`](@ref)- or [`ODD`](@ref)-parity auxiliary density operators (ADOs). The default value of the parameter is `parity=EVEN`.
```julia
Hs::AbstractMatrix  # system Hamiltonian
Bbath::BosonBath    # bosonic   bath object
Fbath::FermionBath  # fermionic bath object
Btier::Int          # bosonic   truncation level 
Ftier::Int          # fermionic truncation level 

# create HEOMLS matrix in EVEN or ODD parity
M_even = M_S(Hs, EVEN)
M_odd  = M_S(Hs, ODD)

M_even = M_Boson(Hs, Btier, Bbath, EVEN) 
M_odd  = M_Boson(Hs, Btier, Bbath, ODD) 

M_even = M_Fermion(Hs, Ftier, Fbath, EVEN) 
M_odd  = M_Fermion(Hs, Ftier, Fbath, ODD)

M_even = M_Boson_Fermion(Hs, Btier, Ftier, Bbath, Fbath, EVEN) 
M_odd  = M_Boson_Fermion(Hs, Btier, Ftier, Bbath, Fbath, ODD) 
```

## Base functions support
### Multiplication between Parity labels
```julia
EVEN * EVEN # gives EVEN
EVEN * ODD  # gives ODD
ODD  * EVEN # gives ODD
ODD  * ODD  # gives EVEN
!EVEN       # gives ODD
!ODD        # gives EVEN
```
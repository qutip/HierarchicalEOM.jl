# [HEOMLS Matrix for Master Equation](@id doc-Master-Equation)

`HEOM.jl` allows users to further add the Lindbladian (superoperator) on any types (`AbstractHEOMMatrix`) of HEOM Liouvillian superoperator ``\hat{\mathcal{M}}``. The Lindbladian describes the dissipative interaction between the system and extra environment.

This method is more efficient than the full HEOM when some of the baths are weakly coupled to the system since it does not require any extra [ADOs](@ref doc-ADOs) space to describe the dynamics and hence reduces the size of ``\hat{\mathcal{M}}``.

## Bosonic Dissipative Env.
If the system is weakly coupled to an extra bosonic environment, the explicit form of the Lindbladian is given by
```math
\hat{\mathcal{D}}(J)\Big[\cdot\Big]=J\left[\cdot\right]J^\dagger-\frac{1}{2}\Big[J^\dagger J, \cdot\Big]_+,
```
where ``J\equiv \sqrt{\gamma}V`` is the jump operator, ``V`` describes the dissipative part (operator) of the dynamics, ``\gamma`` represents a non-negative damping rate and ``[\cdot, \cdot]_+`` stands for anti-commutator.

!!! note "Note"
    The system here can be either bosonic or fermionic. However, if ``V`` is acting on fermionic systems, it should be even-parity to be compatible with charge conservation.

One can add the Lindbladian ``\hat{\mathcal{D}}`` of bosonic environment to the HEOM Liouvillian superoperator ``\hat{\mathcal{M}}`` by calling 

[`addBosonDissipator(M::AbstractHEOMMatrix, jumpOP)`](@ref addBosonDissipator) with the parameters:
 - `M::AbstractHEOMMatrix` : the matrix given from HEOM model
 - `jumpOP::AbstractVector` : The list of collapse (jump) operators ``\{J_i\}_i`` to add. Defaults to empty vector `[]`.

Finally, the function returns a new ``\hat{\mathcal{M}}`` with the same type:
```julia
M0::AbstractHEOMMatrix
J = [J1, J2, ..., Jn]  # jump operators

M1 = addBosonDissipator(M0, J)
```

## Fermionic Dissipative Env.
If the fermionic system is weakly coupled to an extra fermionic environment, the explicit form of the Lindbladian acting on `:even`-parity operators is given by
```math
\hat{\mathcal{D}}_{\textrm{even}}(J)\Big[\cdot\Big]=J\left[\cdot\right]J^\dagger-\frac{1}{2}\Big[J^\dagger J, \cdot\Big]_+,
```
where ``J\equiv \sqrt{\gamma}V`` is the jump operator, ``V`` describes the dissipative part (operator) of the dynamics, ``\gamma`` represents a non-negative damping rate and ``[\cdot, \cdot]_+`` stands for anti-commutator.

For acting on `:odd`-parity operators, the explicit form of the Lindbladian is given by
```math
\hat{\mathcal{D}}_{\textrm{odd}}(J)\Big[\cdot\Big]=-J\left[\cdot\right]J^\dagger-\frac{1}{2}\Big[J^\dagger J, \cdot\Big]_+,
```

One can add the Lindbladian ``\hat{\mathcal{D}}`` of fermionic environment to the HEOM Liouvillian superoperator ``\hat{\mathcal{M}}`` by calling 

[`addFermionDissipator(M::AbstractHEOMMatrix, jumpOP)`](@ref addFermionDissipator) with the parameters:
 - `M::AbstractHEOMMatrix` : the matrix given from HEOM model
 - `jumpOP::AbstractVector` : The list of collapse (jump) operators ``\{J_i\}_i`` to add. Defaults to empty vector `[]`.

!!! note "Parity"
    The parity of the dissipator will be automatically determined by the [parity](@ref doc-Parity) of the given HEOM matrix `M`.

Finally, the function returns a new ``\hat{\mathcal{M}}`` with the same type and parity:
```julia
M0_even::AbstractHEOMMatrix # constructed with :even-parity
M0_odd::AbstractHEOMMatrix  # constructed with :odd-parity
J = [J1, J2, ..., Jn]  # jump operators

M1_even = addFermionDissipator(M0_even, J)
M1_odd  = addFermionDissipator(M0_odd,  J)
```
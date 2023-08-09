# [Auxiliary Density Operators](@id doc-ADOs)

## Introduction
The auxiliary density operators ([`ADOs`](@ref)) ``\rho_{\textbf{j}\vert\textbf{q}}^{(m,n,p)}(t)`` encode environmental effects related to different exponential terms ([`Exponent`](@ref)) present in the [`Bosonic Bath`](@ref doc-Bosonic-Bath) and [`Fermionic Bath`](@ref doc-Fermionic-Bath) correlation functions and provide an iterative description of high-order system-baths memory effects.

In ``\rho_{\textbf{j}\vert\textbf{q}}^{(m,n,p)}(t)``, the tuple ``(m, n, p)`` represents the ``m``th-level-bosonic-and-``n``th-level-fermionic ADO with parity ``p``, and ``\textbf{j}`` (``\textbf{q}``) denotes a vector ``[j_m,\cdots,j_1]`` (``[q_n,\cdots,q_1]``) where each ``j`` (``q``) represents a specific multi-index ensemble ``\{\beta, l\}`` (``\{\alpha, h\}``) with
 - ``\beta`` : denotes the index of bosonic bath
 - ``\alpha`` : denotes the index of fermionic bath
 - ``l`` : denotes the index of exponent in the bosonic bath
 - ``h`` : denotes the index of exponent in the fermionic bath

!!! note "Reduced Density Operator"
    The system reduced density operator refers to ``m=n=0``, namely ``\rho_{\vert}^{(0,0,p)}(t)``.

In `HierarchicalEOM.jl`, we express all the auxiliary density operators into a single column vector and store it in the object defined as : 

[`struct ADOs`](@ref ADOs), 

which is usually obtained after solving the time [evolution](@ref doc-Time-Evolution) or [stationary state](@ref doc-Stationary-State) by a given [HEOM Liouvillian superoperator Matrix](@ref doc-HEOMLS-Matrix).

## Fields
The fields of the structure [`ADOs`](@ref) are as follows:
 - `data` : the vectorized auxiliary density operators
 - `dim` : the dimension of the system
 - `N` : the number of auxiliary density operators
- `parity`: the [parity](@ref doc-Parity) label

One obtain the value of each fields as follows:
```julia
# usually obtained after solving time evolution or stationary state
ados::ADOs

ados.data
ados.dim
ados.N
ados.parity
```
!!! warning "Warning"
    We express all the auxiliary density operators in only a single column vector `ADOs.data`. To obtain each auxiliary density operators in matrix form, please use the following methods and functions.

## Reduced Density Operator
In order to obtain the system reduced density operator in matrix form, just simply call [`getRho`](@ref)
```julia
ados::ADOs
ρ = getRho(ados)
```

## High-Level Auxiliary Density Operators
Although we express all the auxiliary density operators in the vector `ADOs.data`, we still make the [`ADOs`](@ref) like a list where accessing each element would return a specific auxiliary density operator in matrix type. 

In order to obtain the auxiliary density operator in matrix form with a specific index `i`, just simply call [`getADO`](@ref)
```julia
ados::ADOs
i::Int

ρ   = getADO(ados, 1) # the first element will always be the reduced density operator
ado = getADO(ados, i) # the i-th auxiliary density operator
```

Also, [`ADOs`](@ref) supports all the element-wise methods (functions) :

Supports `length(::ADOs)` which returns the total number of auxiliary density operators (same as `ADOs.N`) :
```julia
ados::ADOs
length(ados)
```

Supports bracket operation `[]` which is similar to access the element of a list :
```julia
ados::ADOs

# all the following returned ADO will be in matrix form
ados[1]    # returns the first auxiliary density operator (which is always the reduced density operator)
ados[10]   # returns the 10-th auxiliary density operator
ados[3:10] # returns a list of auxiliary density operators from index 3 to 10
ados[end]  # returns the last auxiliary density operator
```

Supports iteration (`for`-loop) process :
```julia
ados::ADOs

for ado in ados  # iteration
    ado # each auxiliary density operator in matrix form
end
```

!!! note "Work on high-level auxiliary density operators with Hierarchy Dictionary"
    To find the index of the auxiliary density operator and it's corresponding bath [`Exponent`](@ref), please refer to [`Hierarchy Dictionary`](@ref doc-Hierarchy-Dictionary) for more details.

## Expectation Value
Given an observable ``A`` and `ADOs` ``\rho^{(m,n,p)}_{\textbf{j} \vert \textbf{q}}``, one can calculate the expectation value by
```math
\langle A \rangle = \textrm{Tr}\left[A \rho^{(0,0,p)}_{ \vert }\right],
```
where, ``m=n=0`` represents the reduced density operator.

One can directly calculate the expectation values using the function [`Expect`](@ref):
```julia
A::AbstractMatrix # observable

# with a single ADOs
ados::ADOs
E = Expect(A, ados)

# with a list contains many ADOs
ados_list::Vector{ADOs}
Elist = Expect(A, ados_list)
```
Here, `Elist` contains the expectation values corresponding to the `ados_list`.
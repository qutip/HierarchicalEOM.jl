# [Spectrum](@id doc-Spectrum)

## Introduction
We briefly summarize how to numerically compute the spectrum associated with the system degree of freedom. [Phys. Rev. Lett. 109, 266403 (2012)](https://link.aps.org/doi/10.1103/PhysRevLett.109.266403) showed that the spectrum can be evaluated either in time or frequency domain.

`HierarchicalEOM.jl` provides the following listed functions which performs the calculation of spectrum in frequency domain.

 - [Power Spectrum](@ref doc-PS)
 - [Density of States](@ref doc-DOS)

`HierarchicalEOM.jl` wraps some of the functions in [LinearSolve.jl](http://linearsolve.sciml.ai/stable/), which is a very rich numerical library for solving the linear problems and provides many solvers. It offers quite a few options for the user to tailor the solver to their specific needs. The default solver (and its corresponding settings) are chosen to suit commonly encountered problems and should work fine for most of the cases. If you require more specialized methods, such as the choice of algorithm, please refer to [LinearSolve solvers](@ref LS-solvers) and also the documentation of [LinearSolve.jl](http://linearsolve.sciml.ai/stable/).

!!! compat "Extension for CUDA.jl"
    `HierarchicalEOM.jl` provides an extension to support GPU ([`CUDA.jl`](https://github.com/JuliaGPU/CUDA.jl)) acceleration for solving the spectrum, but this feature requires `Julia 1.9+` and `HierarchicalEOM 1.1+`. See [here](@ref doc-ext-CUDA) for more details.

The output of the above listed functions will always be in the type of `Vector{Float64}`, which contains the list of the spectrum values corresponding to the given `ωlist`.

## Common and optional parameters
Furthermore, there are two common optional parameters for all the functions provided below:
- `verbose::Bool` : To display verbose output and progress bar during the process or not. Defaults to `true`.
- `filename::String` : If filename was specified, the value of spectrum for each ω will be saved into the file "filename.txt" during the solving process.

If the filename is specified, the function will automatically save (update) the value (together with a comma behind it) to a new line in the file (with ".txt" behind the filename) once it obtains the solution of each specified ``\omega``. For example, if you specify `filename="test"` and `ωlist=0:1:5`, you will obtain a file `test.txt` where each line in this file (as shown below) is the result of spectrum corresponding to the given `ωlist`:
```
# (the content inside test.txt) #
0.4242990296334028,
0.28617768129333854,
0.21332961856387556,
0.1751179183484055,
0.15739257286986685,
0.1518018484057393,

```

For your convenience, we add those commas (",") in the end of each line for the users to easily do "copy-and-paste" and load these results back into julia's kernel (construct a vector of the given results) again, namely

```julia
results = [
0.4242990296334028,
0.28617768129333854,
0.21332961856387556,
0.1751179183484055,
0.15739257286986685,
0.1518018484057393
]
```

## [Power Spectrum](@id doc-PS)
Start from the power spectrum in the time-domain. We write the system two-time correlation function in terms of the propagator ``\hat{\mathcal{G}}(t)=\exp(\hat{\mathcal{M}} t)`` for ``t>0``. The power spectrum ``\pi S(\omega)`` can be obtained as
```math
\begin{aligned}
\pi S(\omega) 
&= \textrm{Re}\left\{\int_0^\infty dt \langle P(t)Q(0)\rangle e^{-i\omega t}\right\}\\
&= \textrm{Re}\left\{\int_0^\infty dt \langle P e^{\hat{\mathcal{M}} t}Q\rangle e^{-i\omega t}\right\}\\
&= -\textrm{Re}\left\{\langle P (\hat{\mathcal{M}} -i\omega)^{-1} Q\rangle\right\}\\
&= -\textrm{Re}\left\{\textrm{Tr}\left[ P (\hat{\mathcal{M}} -i\omega)^{-1} Q\rho^{(m,n,+)}_{\textbf{j} \vert \textbf{q}}\right]\right\},
\end{aligned}
```
where a half-Fourier transform has been introduced in the third line. We note that only the reduced density operator (``m=n=0``) is considered when taking the final trace operation.

This function solves the linear problem ``\textbf{A x}=\textbf{b}`` at a fixed frequency ``\omega`` where 
 - ``\textbf{A}=\hat{\mathcal{M}}-i\omega``
 - ``\textbf{b}=Q\rho^{(m,n,+)}_{\textbf{j} \vert \textbf{q}}`` 
using the package [LinearSolve.jl](http://linearsolve.sciml.ai/stable/).

Finially, one can obtain the value of the power spectrum for specific ``\omega``, namely
```math
\pi S(\omega) = -\textrm{Re}\left\{\textrm{Tr}\left[ P \textbf{x}\right]\right\}.
```

!!! note "Odd-Parity for Power Spectrum"
    When ``Q`` is an operator acting on fermionic systems and has `ODD`-parity, the HEOMLS matrix ``\hat{\mathcal{M}}`` is acting on the `ODD`-parity space because ``\textbf{b}=Q\rho^{(m,n,+)}_{\textbf{j} \vert \textbf{q}}``. Therefore, remember to construct ``\hat{\mathcal{M}}`` with `ODD` [parity](@ref doc-Parity) in this kind of cases.

See also the docstring : 
```@docs
PowerSpectrum(M::AbstractHEOMLSMatrix, ρ::Union{QuantumObject,ADOs}, P_op::Union{QuantumObject,HEOMSuperOp}, Q_op::Union{QuantumObject,HEOMSuperOp}, ωlist::AbstractVector, reverse::Bool = false; solver = UMFPACKFactorization(), verbose::Bool = true, filename::String = "", SOLVEROptions...)
```

```julia
M::AbstractHEOMLSMatrix

# the input state can be in either type (but usually ADOs):
ρ::QuantumObject # the reduced density operator
ρ::ADOs # the ADOs solved from "evolution" or "steadystate"

P::QuantumObject 
Q::QuantumObject

# the spectrum value for the specific frequency ω which need to be solved
ωlist = 0:0.5:2 # [0.0, 0.5, 1.0, 1.5, 2.0]

πSω = spectrum(M, ρ,    Q, ωlist) # P will automatically be considered as "the adjoint of Q operator"
πSω = spectrum(M, ρ, P, Q, ωlist) # user specify both P and Q operator
```

## [Density of States](@id doc-DOS)
Start from the density of states for fermionic systems in the time-domain. We write the system two-time correlation function in terms of the propagator ``\hat{\mathcal{G}}(t)=\exp(\hat{\mathcal{M}} t)`` for ``t>0``. The density of states ``\pi A(\omega)`` can be obtained as
```math
\begin{aligned}
\pi A(\omega) 
&= \textrm{Re}\left\{\int_0^\infty dt \langle d(t)d^\dagger(0)\rangle e^{i\omega t}\right\} + \textrm{Re}\left\{\int_0^\infty dt \langle d^\dagger(t)d(0)\rangle e^{-i\omega t}\right\}\\
&= \textrm{Re}\left\{\int_0^\infty dt \langle d e^{\hat{\mathcal{M}} t}d^\dagger\rangle e^{i\omega t}\right\}+\textrm{Re}\left\{\int_0^\infty dt \langle d^\dagger e^{\hat{\mathcal{M}} t}d\rangle e^{-i\omega t}\right\}\\
&= -\textrm{Re}\left\{\langle d (\hat{\mathcal{M}} +i\omega)^{-1} d^\dagger\rangle + \langle d^\dagger (\hat{\mathcal{M}} -i\omega)^{-1} d\rangle\right\}\\
&= -\textrm{Re}\left\{\textrm{Tr}\left[ d (\hat{\mathcal{M}} +i\omega)^{-1} d^\dagger\rho^{(m,n,+)}_{\textbf{j} \vert \textbf{q}}\right] + \textrm{Tr}\left[ d^\dagger (\hat{\mathcal{M}} -i\omega)^{-1} d\rho^{(m,n,+)}_{\textbf{j} \vert \textbf{q}}\right]\right\},
\end{aligned}
```
where a half-Fourier transform has been introduced in the third line. We note that only the reduced density operator (``m=n=0``) is considered when taking the final trace operation.

This functionsolves two linear problems ``\textbf{A}_+ \textbf{x}_+=\textbf{b}_+`` and ``\textbf{A}_- \textbf{x}_-=\textbf{b}_-`` at a fixed frequency ``\omega`` where 
 - ``\textbf{A}_+=\hat{\mathcal{M}}+i\omega``
 - ``\textbf{b}_+=d^\dagger\rho^{(m,n,+)}_{\textbf{j} \vert \textbf{q}}`` 
 - ``\textbf{A}_-=\hat{\mathcal{M}}-i\omega``
 - ``\textbf{b}_-=d\rho^{(m,n,+)}_{\textbf{j} \vert \textbf{q}}`` 
using the package [LinearSolve.jl](http://linearsolve.sciml.ai/stable/).

Finially, one can obtain the density of states for specific ``\omega``, namely
```math
\pi A(\omega) = -\textrm{Re}\left\{\textrm{Tr}\left[ d \textbf{x}_+\right]+\textrm{Tr}\left[ d^\dagger \textbf{x}_-\right]\right\}.
```

!!! note "Odd-Parity for Density of States"
    As shown above, the HEOMLS matrix ``\hat{\mathcal{M}}`` acts on the `ODD`-parity space, compatibly with the parity of both the operators ``\textbf{b}_-=d\rho^{(m,n,+)}_{\textbf{j} \vert \textbf{q}}`` and ``\textbf{b}_+=d^\dagger\rho^{(m,n,+)}_{\textbf{j} \vert \textbf{q}}``. Therefore, remember to construct ``\hat{\mathcal{M}}`` with `ODD` [parity](@ref doc-Parity) for solving spectrum of fermionic systems.

See also the docstring : 
```@docs
DensityOfStates(M::AbstractHEOMLSMatrix, ρ::QuantumObject, d_op::QuantumObject, ωlist::AbstractVector; solver=UMFPACKFactorization(), verbose::Bool = true, filename::String = "", SOLVEROptions...)
```

```julia
Hs::QuantumObject  # system Hamiltonian
bath::FermionBath  # fermionic bath object
tier::Int          # fermionic truncation level 

# create HEOMLS matrix in both :even and ODD parity
M_even = M_Fermion(Hs, tier, bath) 
M_odd  = M_Fermion(Hs, tier, bath, ODD) 

# the input state can be in either type of density operator matrix or ADOs (but usually ADOs):
ados = steadystate(M_even)

# the (usually annihilation) operator "d" as shown above
d::QuantumObject 

# the spectrum value for the specific frequency ω which need to be solved
ω_list = 0:0.5:2 # [0.0, 0.5, 1.0, 1.5, 2.0]

πAω = DensityOfStates(M_odd, ados, d, ω_list)
```
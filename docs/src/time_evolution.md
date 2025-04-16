# [Time Evolution](@id doc-Time-Evolution)
## Introduction
`HierarchicalEOM.jl` implements various methods and solvers to simulate the open quantum system dynamics. 
The [HEOM Liouvillian superoperator (HEOMLS) matrix](@ref doc-HEOMLS-Matrix) ``\hat{\mathcal{M}}`` characterizes the dynamics of the reduce state and in the full extended space of all [auxiliary density operators (ADOs)](@ref doc-ADOs) ``\rho^{(m,n,p)}_{\textbf{j} \vert \textbf{q}}(t)``, namely
```math
\begin{equation}
\partial_{t}\rho^{(m,n,p)}_{\textbf{j} \vert \textbf{q}}(t)=\hat{\mathcal{M}}\rho^{(m,n,p)}_{\textbf{j} \vert \textbf{q}}(t)
\end{equation}
```

### `HEOMsolve` and `TimeEvolutionHEOMSol`
To solve the dynamics of the reduced state and also all the [ADOs](@ref doc-ADOs), you only need to call [`HEOMsolve`](@ref). Different methods (see the contents below) are implemented with different input parameters of the function which makes it easy to switch between different methods. The output of the function [`HEOMsolve`](@ref) for each methods will always be in the type [`TimeEvolutionHEOMSol`](@ref), which contains the results (including [`ADOs`](@ref) and expectation values at each time point) and some information from the solver. One can obtain the value of each fields in [`TimeEvolutionHEOMSol`](@ref) as follows:

```julia
sol::TimeEvolutionHEOMSol

sol.Btier   # the tier (cutoff level) for bosonic hierarchy
sol.Ftier   # the tier (cutoff level) for fermionic hierarchy
sol.times   # The time list of the evolution.
sol.ados    # The list of result ADOs at each time point.
sol.expect  # The expectation values corresponding to each time point in `times`.
sol.retcode # The return code from the solver.
sol.alg     # The algorithm which is used during the solving process.
sol.abstol  # The absolute tolerance which is used during the solving process.
sol.reltol  # The relative tolerance which is used during the solving process.
```

### Expectation Values
Given an observable ``A`` and the [`ADOs`](@ref doc-ADOs) ``\rho^{(m,n,p)}_{\textbf{j} \vert \textbf{q}}(t)``, one can calculate the expectation value by
```math
\langle A(t) \rangle = \textrm{Tr}\left[A \rho^{(0,0,p)}_{ \vert }(t)\right],
```
where, ``m=n=0`` represents the reduced density operator, see [`ADOs`](@ref doc-ADOs) for more details.

One can directly calculate the expectation value by specifying the keyword argument `e_ops` (a list of observables), and the expectation values corresponding to each time point and observables will be stored in [`TimeEvolutionHEOMSol`](@ref):
```julia
A1::QuantumObject # observable 1
A2::QuantumObject # observable 2
sol = HEOMsolve(...; e_ops = [A1, A2], ...) # the input parameters depend on the different methods you choose.
sol.expect[1,:] # the expectation values of observable 1 (`A1`) corresponding to each time point in `sol.times`
sol.expect[2,:] # the expectation values of observable 2 (`A2`) corresponding to each time point in `sol.times`
```

An alternative way for calculating the expectation values is to use the function [`QuantumToolbox.expect`](@ref) together with the list of [`ADOs`](@ref) stored in [`TimeEvolutionHEOMSol`](@ref):
```julia
A::QuantumObject # observable
sol = HEOMsolve(...) # the input parameters depend on the different methods you choose.
ados_list = sol.ados

Elist = expect(A, ados_list)
```
Here, `Elist` contains the expectation values corresponding to the `ados_list` (i.e., the reduced density operator in each time step).

## Common and optional parameters
There are three common optional parameters for all the methods provided below:
 - `e_ops::Union{Nothing,AbstractVector}`: List of operators for which to calculate expectation values.
 - `verbose::Bool` : To display verbose output and progress bar during the process or not. Defaults to `true`.
 - `filename::String` : If filename was specified, the ADOs at each time point will be saved into the JLD2 file after the solving process. Default to Empty `String`: `""`.

If the `filename` is specified, the function will automatically save the [`ADOs`](@ref) to the file (with `.jld2` behind the `filename`) once the solving process is finished. The saving method is based on the package [`JLD2.jl`](https://juliaio.github.io/JLD2.jl/stable/), which saves and loads `Julia` data structures in a format comprising a subset of HDF5.
```julia
tlist = 0:0.5:5
ados_list = HEOMsolve(..., tlist, ...; filename="test", ...)
```
The solution of the [`ADOs`](@ref) for each time step in `tlist` is saved in the file named `test.jld2` with a key: `"ados"`.

To retrieve the solution the list of [`ADOs`](@ref) from a previously saved file `"text.jld2"`, just read the file with the methods provided by [`JLD2.jl`](https://juliaio.github.io/JLD2.jl/stable/) and specify the key: `"ados"`, namely
```julia
using HierarchicalEOM, JLD2 # remember to import these before retrieving the solution

filename = "test.jld2"
jldopen(filename, "r") do file
    ados_list = file["ados"]
end
```

## Ordinary Differential Equation (ODE) Method
The first method is implemented by solving the ordinary differential equation (ODE). `HierarchicalEOM.jl` wraps some of the functions in [`DifferentialEquations.jl`](https://diffeq.sciml.ai/stable/), which is a very rich numerical library for solving the differential equations and provides many ODE solvers. It offers quite a few options for the user to tailor the solver to their specific needs. The default solver (and its corresponding settings) are chosen to suit commonly encountered problems and should work fine for most of the cases. If you require more specialized methods, such as the choice of algorithm, please refer to [DifferentialEquations solvers](@ref ODE-solvers) and also the documentation of [`DifferentialEquations.jl`](https://diffeq.sciml.ai/stable/).

!!! compat "Extension for CUDA.jl"
    `HierarchicalEOM.jl` provides an extension to support GPU ([`CUDA.jl`](https://github.com/JuliaGPU/CUDA.jl)) acceleration for [`HEOMsolve`](@ref) (only for ODE method). See [here](@ref doc-ext-CUDA) for more details.

See the docstring of this method:  

```@docs
HEOMsolve(M::AbstractHEOMLSMatrix, ρ0::T_state, tlist::AbstractVector; e_ops::Union{Nothing,AbstractVector} = nothing, solver::OrdinaryDiffEqAlgorithm = DP5(), H_t::Union{Nothing,Function} = nothing, params::NamedTuple = NamedTuple(), verbose::Bool = true, filename::String = "", SOLVEROptions...,) where {T_state<:Union{QuantumObject,ADOs}}
```

```julia
# the time-independent HEOMLS matrix
M::AbstractHEOMLSMatrix  

# the initial state can be either the system density operator or ADOs
ρ0::QuantumObject
ρ0::ADOs

# specific time points to save the solution during the solving process.  
tlist = 0:0.5:2 # [0.0, 0.5, 1.0, 1.5, 2.0]

sol = HEOMsolve(M, ρ0, tlist)
```

## Time Dependent Problems
In general, the time-dependent system Hamiltonian can be separated into the time-independent and time-dependent parts, namely
```math
H_s (t) = H_0 + H_1(t).
```
We again wrap some of the functions in [`DifferentialEquations.jl`](https://diffeq.sciml.ai/stable/) to solve the time-dependent problems here.

To deal with the time-dependent system Hamiltonian problem in `HierarchicalEOM.jl`, we first construct the [HEOMLS matrices](@ref doc-HEOMLS-Matrix) ``\hat{\mathcal{M}}`` with **time-independent** Hamiltonian ``H_0``:

```julia
M = M_S(H0, ...)
M = M_Boson(H0, ...)
M = M_Fermion(H0, ...)
M = M_BosonFermion(H0, ...)
```

To solve the dynamics characterized by ``\hat{\mathcal{M}}`` together with the time-dependent part of system Hamiltonian ``H_1(t)``, you can specify keyword arguments `H_t` and `params` while calling [`HEOMsolve`](@ref). Here, `H_t` must be specified as a `QuantumToolbox.QuantumObjectEvolution` (or `QobjEvo`), and `params` should contain all the extra parameters you need for `QobjEvo`, for example:
```julia
# in this case, p will be passed in as a NamedTuple: (p0 = p0, p1 = p1, p2 = p2)
coef(p::NamedTuple, t) = sin(p.p0 * t) + sin(p.p1 * t) + sin(p.p2 * t)

σx = sigmax() # Pauli-X matrix
H_1 = QobjEvo(σx, coef)
```
The `p` can be passed to `H_1` directly from the keyword argument in [`HEOMsolve`](@ref) called `params`:
```julia
M::AbstractHEOMLSMatrix
ρ0::QuantumObject
tlist = 0:0.1:10
p = (p0 = 0.1, p1 = 1, p2 = 10)

sol = HEOMsolve(M, ρ0, tlist; H_t = H_1, params = p)
```

## Propagator Method
The second method is implemented by directly construct the propagator of a given [HEOMLS matrix](@ref doc-HEOMLS-Matrix) ``\hat{\mathcal{M}}``. Because ``\hat{\mathcal{M}}`` is time-independent, the equation above can be solved analytically as
```math
\rho^{(m,n,p)}_{\textbf{j} \vert \textbf{q}}(t)=\hat{\mathcal{G}}(t)\rho^{(m,n,p)}_{\textbf{j} \vert \textbf{q}}(0),
```
where ``\hat{\mathcal{G}}(t)\equiv \exp(\hat{\mathcal{M}}t)`` is the propagator for all [ADOs](@ref doc-ADOs) corresponding to ``\hat{\mathcal{M}}``.

To construct the propagator, we wrap the function in the package [`fastExpm.jl`](https://github.com/fmentink/FastExpm.jl), which is optimized for the exponentiation of either large-dense or sparse matrices.

See the docstring of this method:  

```@docs
HEOMsolve(M::AbstractHEOMLSMatrix ,ρ0::T_state, Δt::Real, steps::Int; e_ops::Union{Nothing,AbstractVector} = nothing, threshold = 1.0e-6, nonzero_tol = 1.0e-14, verbose::Bool = true, filename::String = "",) where {T_state<:Union{QuantumObject,ADOs}}
```

```julia
# the time-independent HEOMLS matrix
M::AbstractHEOMLSMatrix  

# the initial state can be either the system density operator or ADOs
ρ0::QuantumObject
ρ0::ADOs

# A specific time interval (time step)
Δt = 0.5

# The number of time steps for the propagator to apply
steps = 4

# equivalent to tlist = 0 : Δt : (Δt * steps)
sol = HEOMsolve(M, ρ0, Δt, steps) 
```

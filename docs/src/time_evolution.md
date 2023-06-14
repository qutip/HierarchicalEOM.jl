# [Time Evolution](@id doc-Time-Evolution)
## Introduction
`HEOM.jl` implements various methods and solvers to simulate the open quantum system dynamics. 
The [HEOM Liouvillian superoperator (HEOMLS) matrix](@ref doc-HEOMLS-Matrix) ``\hat{\mathcal{M}}`` characterizes the dynamics of the reduce state and in the full extended space of all [auxiliary density operators (ADOs)](@ref doc-ADOs) ``\rho^{(m,n,p)}_{\textbf{j} \vert \textbf{q}}(t)``, namely
```math
\begin{equation}
\partial_{t}\rho^{(m,n,p)}_{\textbf{j} \vert \textbf{q}}(t)=\hat{\mathcal{M}}\rho^{(m,n,p)}_{\textbf{j} \vert \textbf{q}}(t)
\end{equation}
```

### Output of `evolution`
To solve the dynamics of the reduced state and also all the [ADOs](@ref doc-ADOs), you only need to call [`evolution`](@ref). Different methods are implemented with different input parameters of the function which makes it easy to switch between different methods. The output of the function [`evolution`](@ref) for each methods will always be in the type `Vector{ADOs}`, which contains a list of [`ADOs`](@ref) corresponds to the given time steps.

### Expectation Values
Given an observable ``A`` and the [`ADOs`](@ref doc-ADOs) ``\rho^{(m,n,p)}_{\textbf{j} \vert \textbf{q}}(t)``, one can calculate the expectation value by
```math
\langle A(t) \rangle = \textrm{Tr}\left[A \rho^{(0,0,p)}_{ \vert }(t)\right],
```
where, ``m=n=0`` represents the reduced density operator, see [`ADOs`](@ref doc-ADOs) for more details.

One can directly calculate the expectation values using the function [`Expect`](@ref) together with the output of [`evolution`](@ref):
```julia
A::AbstractMatrix # observable
ados_list = evolution(...) # the input parameters depend on the different methods you choose.

Elist = Expect(A, ados_list)
```
Here, `Elist` contains the expectation values corresponding to the `ados_list` (i.e., the reduced density operator in each time step).

### Common and optional parameters for `evolution`
Furthermore, there are two common optional parameters for all the methods provided below:
 - `verbose::Bool` : To display verbose output and progress bar during the process or not. Defaults to `true`.
 - `filename::String` : If filename was specified, the ADOs at each time point will be saved into the JLD2 file during the solving process. Default to Empty `String`: `""`.
If the filename is specified, the function will automatically save (update) the [`ADOs`](@ref) to the file (with ".jld2" behind the filename) once it obtains the solution of the specified time step. The saving method is based on the package [`JLD2.jl`](https://juliaio.github.io/JLD2.jl/stable/), which saves and loads `Julia` data structures in a format comprising a subset of HDF5.
```julia
tlist = 0:0.5:5
ados_list = evolution(..., tlist, ...; filename="test", ...)
```
The solution of the ADOs for each time step in `tlist` is saved in the file named `"test.jld2"`.

To retrieve the solution ([`ADOs`](@ref)) from a previously saved file `"text.jld2"`, just read the file with the methods provided by [`JLD2.jl`](https://juliaio.github.io/JLD2.jl/stable/). The solution for a specific time step can be extract by using the string of the time step as the `"key"`.
For example, if you want to obtain the solution at time `1.5`, which is one of the time steps in `tlist`:
```julia
using HEOM, JLD2 # rembember to import these before retrieving the solution

t = 1.5
filename = "test.jld2"
jldopen(filename, "r") do file
    ados = file[string(t)]
end
```

## Ordinary Differential Equation Method
The first method is implemented by solving the ordinary differential equation (ODE) as shown above. `HEOM.jl` wraps some of the functions in [`DifferentialEquations.jl`](https://diffeq.sciml.ai/stable/), which is a very rich numerical library for solving the differential equations and provides many ODE solvers. It offers quite a few options for the user to tailor the solver to their specific needs. The default solver (and its corresponding settings) are chosen to suit commonly encountered problems and should work fine for most of the cases. If you require more specialized methods, such as the choice of algorithm, please refer to the documentation of [`DifferentialEquations.jl`](https://diffeq.sciml.ai/stable/).

### Given the initial state as Density Operator (`AbstractMatrix` type)

See the docstring of this method:  

[`evolution(M::AbstractHEOMMatrix, ρ0, tlist::AbstractVector)`](@ref)

```julia
# the time-independent HEOMLS matrix
M::AbstractHEOMMatrix  

# the initial state of the system density operator
ρ0::AbstractMatrix

# specific time points to save the solution during the solving process.  
tlist = 0:0.5:2 # [0.0, 0.5, 1.0, 1.5, 2.0]

ados_list = evolution(M, ρ0, tlist)
```

### Given the initial state as Auxiliary Density Operators
!!! note "Note" 
    This method is usually used when you want to solve the time evolution again with the initial state are given from the last time point of the previous result.
See the docstring of this method:   

[`evolution(M::AbstractHEOMMatrix, ados::ADOs, tlist::AbstractVector)`](@ref)

```julia
# the time-independent HEOMLS matrix
M::AbstractHEOMMatrix  

# the initial state of the ADOs (usually obtianed from previous solving result)
ados::ADOs      

# specific time points to save the solution during the solving process. 
tlist = 0:0.5:2 # [0.0, 0.5, 1.0, 1.5, 2.0]

ados_list = evolution(M, ados, tlist)
```

## Propagator Method
The second method is implemented by directly construct the propagator of a given [HEOMLS matrix](@ref doc-HEOMLS-Matrix) ``\hat{\mathcal{M}}``. Because ``\hat{\mathcal{M}}`` is time-independent, the equation above can be solved analytically as
```math
\rho^{(m,n,p)}_{\textbf{j} \vert \textbf{q}}(t)=\hat{\mathcal{G}}(t)\rho^{(m,n,p)}_{\textbf{j} \vert \textbf{q}}(0),
```
where ``\hat{\mathcal{G}}(t)\equiv \exp(\hat{\mathcal{M}}t)`` is the propagator for all [ADOs](@ref doc-ADOs) corresponding to ``\hat{\mathcal{M}}``.

To construct the propagator, we wrap the function in the package [`fastExpm.jl`](https://github.com/fmentink/FastExpm.jl), which is optimized for the exponentiation of either large-dense or sparse matrices.

### Given the initial state as Density Operator (`AbstractMatrix` type)
See the docstring of this method:  

[`evolution(M::AbstractHEOMMatrix, ρ0, Δt::Real, steps::Int)`](@ref)

```julia
# the time-independent HEOMLS matrix
M::AbstractHEOMMatrix  

# the initial state of the system density operator
ρ0::AbstractMatrix

# A specific time interval (time step)
Δt = 0.5

# The number of time steps for the propagator to apply
steps = 4

# equivalent to tlist = 0 : Δt : (Δt * steps)
ados_list = evolution(M, ρ0, Δt, steps) 
```

### Given the initial state as Auxiliary Density Operators
!!! note "Note" 
    This method is usually used when you want to solve the time evolution again with the initial state are given from the last time point of the previous result.
See the docstring of this method:  

[`evolution(M::AbstractHEOMMatrix, ados::ADOs, Δt::Real, steps::Int)`](@ref)

```julia
# the time-independent HEOMLS matrix
M::AbstractHEOMMatrix  

# the initial state of the ADOs (usually obtianed from previous solving result)
ados::ADOs

# A specific time interval (time step)
Δt = 0.5

# The number of time steps for the propagator to apply
steps = 4

# equivalent to tlist = 0 : Δt : (Δt * steps)
ados_list = evolution(M, ados, Δt, steps) 
```

## Time Dependent Problems
In general, the time-dependent system Hamiltonian can be separated into the time-independent and time-dependent parts, namely
```math
H_s (t) = H_0 + H_1(t).
```
We again wrap some of the functions in [`DifferentialEquations.jl`](https://diffeq.sciml.ai/stable/) to solve the time-dependent problems here.

To deal with the time-dependent system Hamiltonian problem in `HEOM.jl`, we first construct the [HEOMLS matrices](@ref doc-HEOMLS-Matrix) ``\hat{\mathcal{M}}`` with **time-independent** Hamiltonian ``H_0``:
```julia
M = M_S(H0, ...)
M = M_Boson(H0, ...)
M = M_Fermion(H0, ...)
M = M_BosonFermion(H0, ...)
```
To solve the dynamics characterized by ``\hat{\mathcal{M}}`` together with the time-dependent part of system Hamiltonian ``H_1(t)``, you can call either of the following two functions (one takes the type of initial state as density matrix and the other one takes [`ADOs`](@ref)):
 - [`evolution(M::AbstractHEOMMatrix, ρ0, tlist::AbstractVector, H::Function, param::Tuple = ())`](@ref)
 - [`evolution(M::AbstractHEOMMatrix, ados::ADOs, tlist::AbstractVector, H::Function, param::Tuple = ())`](@ref).

Here, the definition of user-defined function `H` must be in the form `H(p::Tuple, t)` and returns the time-dependent part of system Hamiltonian (in `AbstractMatrix` type) at any given time point `t`. The parameter `p` should be a `Tuple` which contains all the extra parameters you need for the function `H`. For example:
```julia
function H_pump(p, t)  
    p0, p1, p2 = p 
    # in this case, p should be passed in as a tuple: (p0, p1, p2) 
    
    σx = [0 1; 1 0] # Pauli-X matrix
    return (sin(p0 * t) + sin(p1 * t) + sin(p2 * t)) * σx
end
```
The parameter tuple `p` will be passed to your function `H` directly from one of the parameter in `evolution` called `param`:
```julia
M::AbstractHEOMMatrix
ρ0::AbstractMatrix
tlist = 0:0.1:10
p = (0.1, 1, 10)

ados_list = evolution(M, ρ0, tlist, H_pump, p)
```

!!! warning "Warning"
    If you don't need any extra `param` in your case, you still need to put a redundant one in the definition of `H`, for example:

```julia
function H_pump(p, t)
    σx = [0 1; 1 0] # Pauli-X matrix
    return sin(0.1 * t) * σx
end

M::AbstractHEOMMatrix
ρ0::AbstractMatrix
tlist = 0:0.1:10

ados_list = evolution(M, ρ0, tlist, H_pump)
```
!!! note "Note"
    The default value for `param` in `evolution` is an empty tuple `()`.

### Given the initial state as Density Operator (`AbstractMatrix` type)
See the docstring of this method:  

[`evolution(M::AbstractHEOMMatrix, ρ0, tlist::AbstractVector, H::Function, param::Tuple = ())`](@ref)


### Given the initial state as Auxiliary Density Operators
!!! note "Note" 
    This method is usually used when you want to solve the time evolution again with the initial state are given from the last time point of the previous result.

See the docstring of this method:  
[`evolution(M::AbstractHEOMMatrix, ados::ADOs, tlist::AbstractVector, H::Function, param::Tuple = ())`](@ref)
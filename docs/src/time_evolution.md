# [Time Evolution](@id doc-Time-Evolution)
`Heom.jl` implements various methods and solvers to simulate the open quantum system dynamics. 
The [HEOM Liouvillian superoperator (HEOMLS) matrix](@ref doc-HEOMLS-Matrix) ``\hat{\mathcal{M}}`` characterizes the dynamics of the reduce state and in the full extended space of all [auxiliary density operators (ADOs)](@ref doc-ADOs) ``\rho^{(m,n,p)}_{\textbf{j} \vert \textbf{q}}(t)``, namely
```math
\begin{equation}
\partial_{t}\rho^{(m,n,p)}_{\textbf{j} \vert \textbf{q}}(t)=\hat{\mathcal{M}}\rho^{(m,n,p)}_{\textbf{j} \vert \textbf{q}}(t)
\end{equation}
```

To solve the dynamics of the reduced state and also all the [ADOs](@ref doc-ADOs), you only need to call [`evolution`](@ref). Different methods are implemented with different input parameters of the function which makes it easy to switch between different methods. The output of the function [`evolution`](@ref) for each methods will always be in the type `Vector{ADOs}`, which contains a list of [`ADOs`](@ref) corresponds to the given time steps.

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
using Heom, JLD2 # rembember to import these before retrieving the solution

t = 1.5
filename = "test.jld2"
jldopen(filename, "r") do file
    ados = file[string(t)]
end
```

## Ordinary Differential Equation Method
The first method is implemented by solving the ordinary differential equation (ODE) as shown above. `Heom.jl` wraps some of the functions in [`DifferentialEquations.jl`](https://diffeq.sciml.ai/stable/), which is a very rich numerical library for solving the differential equations and provides many ODE solvers. It offers quite a few options for the user to tailor the solver to their specific needs. The default solver (and its corresponding settings) are chosen to suit commonly encountered problems and should work fine for most of the cases. If you require more specialized methods, such as the choice of algorithm, please refer to the documentation of [`DifferentialEquations.jl`](https://diffeq.sciml.ai/stable/).

### Given the initial state as Density Operator (`AbstractMatrix` type)

See the docstring of this method:  

[`evolution(M::AbstractHEOMMatrix, ρ0, tlist::AbstractVector; solver = DP5(), reltol::Real = 1.0e-6, abstol::Real = 1.0e-8, maxiters::Real = 1e5, save_everystep::Bool=false, verbose::Bool = true, filename::String = "", SOLVEROptions...)`](@ref)

```julia
# the time-independent HEOMLS matrix
M::AbstractHEOMMatrix  

# the initial state of the system density operator
ρ0::AbstractMatrix

# specific time points to save the solution during the solving process.  
tlist::AbstractVector  

ados_list = (M, ρ0, tlist)
```

### Given the initial state as Auxiliary Density Operators
This method is usually used when you want to solve the time evolution again with the initial state are given from the last time point of the previous result.
See the docstring of this method:   

[`evolution(M::AbstractHEOMMatrix, ados::ADOs, tlist::AbstractVector; solver = DP5(), reltol::Real = 1.0e-6, abstol::Real = 1.0e-8, maxiters::Real = 1e5, save_everystep::Bool=false, verbose::Bool = true, filename::String = "", SOLVEROptions...)`](@ref)

```julia
# the time-independent HEOMLS matrix
M::AbstractHEOMMatrix  

# the initial state of the ADOs (usually obtianed from previous solving result)
ados::ADOs      

# specific time points to save the solution during the solving process. 
tlist::AbstractVector  

ados_list = (M, ados, tlist)
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

[`evolution(M::AbstractHEOMMatrix, ρ0, Δt::Real, steps::Int; threshold   = 1.0e-6, nonzero_tol = 1.0e-14, verbose::Bool = true, filename::String = "")`](@ref)

```julia
# the time-independent HEOMLS matrix
M::AbstractHEOMMatrix  

# the initial state of the system density operator
ρ0::AbstractMatrix

# A specific time interval (time step)
Δt::Real

# The number of time steps for the propagator to apply
steps::Int

ados_list = (M, ρ0, Δt, steps)
```

### Given the initial state as Auxiliary Density Operators
This method is usually used when you want to solve the time evolution again with the initial state are given from the last time point of the previous result.
See the docstring of this method:  

[`evolution(M::AbstractHEOMMatrix, ados::ADOs, Δt::Real, steps::Int; threshold   = 1.0e-6, nonzero_tol = 1.0e-14, verbose::Bool = true, filename::String = "")`](@ref)

```julia
# the time-independent HEOMLS matrix
M::AbstractHEOMMatrix  

# the initial state of the ADOs (usually obtianed from previous solving result)
ados::ADOs

# A specific time interval (time step)
Δt::Real

# The number of time steps for the propagator to apply
steps::Int

ados_list = (M, ρ0, Δt, steps)
```

## Time Dependent Problems

### Given the initial state as Density Operator (`AbstractMatrix` type)
See the docstring of this method:  

[`evolution(M::AbstractHEOMMatrix, ρ0, tlist::AbstractVector, H::Function, param::Tuple = (); solver = DP5(), reltol::Real = 1.0e-6, abstol::Real = 1.0e-8, maxiters::Real = 1e5, save_everystep::Bool=false, verbose::Bool = true, filename::String = "", SOLVEROptions...)`](@ref)


### Given the initial state as Auxiliary Density Operators
This method is usually used when you want to solve the time evolution again with the initial state are given from the last time point of the previous result
See the docstring of this method:  

[`evolution(M::AbstractHEOMMatrix, ados::ADOs, tlist::AbstractVector, H::Function, param::Tuple = (); solver = DP5(), reltol::Real = 1.0e-6, abstol::Real = 1.0e-8, maxiters::Real = 1e5, save_everystep::Bool=false, verbose::Bool = true, filename::String = "", SOLVEROptions...)`](@ref)
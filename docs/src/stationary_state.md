# [Stationary State](@id doc-Stationary-State)
`HierarchicalEOM.jl` implements two different ways to calculate stationary states of all [Auxiliary Density Operators (ADOs)](@ref doc-ADOs).

To solve the stationary state of the reduced state and also all the [ADOs](@ref doc-ADOs), you only need to call [`steadystate`](@ref). Different methods are implemented with different input parameters of the function which makes it easy to switch between different methods. The output of the function [`steadystate`](@ref) for each methods will always be in the type of the auxiliary density operators [`ADOs`](@ref).

## Solve with [LinearSolve.jl](http://linearsolve.sciml.ai/stable/)
The first method is implemented by solving the linear problem
```math
0=\hat{\mathcal{M}}\rho^{(m,n,p)}_{\textbf{j} \vert \textbf{q}}(t)
```

`HierarchicalEOM.jl` wraps some of the functions in [LinearSolve.jl](http://linearsolve.sciml.ai/stable/), which is a very rich numerical library for solving the linear problems and provides many solvers. It offers quite a few options for the user to tailor the solver to their specific needs. The default solver (and its corresponding settings) are chosen to suit commonly encountered problems and should work fine for most of the cases. If you require more specialized methods, such as the choice of algorithm, please refer to [LinearSolve solvers](@ref LS-solvers) and also the documentation of [LinearSolve.jl](http://linearsolve.sciml.ai/stable/).

```julia
# the HEOMLS matrix
M::AbstractHEOMLSMatrix  
ados_steady = steadystate(M)
```
!!! warning "Unphysical solution"
    This method does not require an initial condition ``\rho^{(m,n,p)}_{\textbf{j} \vert \textbf{q}}(0)``. Although this method works for most of the cases, it does not guarantee that one can obtain a physical (or unique) solution. If there is any problem within the solution, please try the second method which solves with an initial condition, as shown below.

## Solve with [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/)
The second method is implemented by solving the ordinary differential equation (ODE) method :
```math
\partial_{t}\rho^{(m,n,p)}_{\textbf{j} \vert \textbf{q}}(t)=\hat{\mathcal{M}}\rho^{(m,n,p)}_{\textbf{j} \vert \textbf{q}}(t)
```
until finding a stationary solution.

`HierarchicalEOM.jl` wraps some of the functions in [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/), which is a very rich numerical library for solving the differential equations and provides many ODE solvers. It offers quite a few options for the user to tailor the solver to their specific needs. The default solver (and its corresponding settings) are chosen to suit commonly encountered problems and should work fine for most of the cases. If you require more specialized methods, such as the choice of algorithm, please refer to the documentation of [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/).

### Given the initial state as Density Operator (`QuantumObject` type)

```julia
# the HEOMLS matrix
M::AbstractHEOMLSMatrix  

# the initial state of the system density operator
ρ0::QuantumObject

ados_steady = steadystate(M, ρ0)
```

### Given the initial state as Auxiliary Density Operators

```julia
# the HEOMLS matrix
M::AbstractHEOMLSMatrix  

# the initial state of the ADOs
ados::ADOs

ados_steady = steadystate(M, ados)
```
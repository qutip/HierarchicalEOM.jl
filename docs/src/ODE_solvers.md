# [DifferentialEquations solvers](@id ODE-solvers)

In this page, we list several recommended `alg`orithms provided by [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/) for solving time evolution in hierarchical equations of motion approach.  

Remember to import `OrdinaryDiffEq.jl` (or `DifferentialEquations.jl`)

```julia
using OrdinaryDiffEq ## or "using DifferentialEquations" 
```

(click [here](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/) to see the full `alg`orithm list provided by `DifferentialEquations.jl`)

For any extra solver options, we can add it in the function `HEOMsolve` with keyword arguments. These keyword arguments will be directly pass to the solvers in `DifferentialEquations`
(click [here](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/) to see the documentation for the common solver options)

### DP5 (Default algorithm)
Dormand-Prince's 5/4 Runge-Kutta method. (free 4th order interpolant)

```julia
DP5()
```

### RK4
The canonical Runge-Kutta Order 4 method. Uses a defect control for adaptive stepping using maximum error over the whole interval.

```julia
RK4()
```

### Tsit5
Tsitouras 5/4 Runge-Kutta method. (free 4th order interpolant).

```julia
Tsit5()
```

### Vern7
Verner's “Most Efficient” 7/6 Runge-Kutta method. (lazy 7th order interpolant).

```julia
Vern7()
```

### Vern9
Verner's “Most Efficient” 9/8 Runge-Kutta method. (lazy 9th order interpolant)

```julia
Vern9()
```

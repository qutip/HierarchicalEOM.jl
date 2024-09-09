# [DifferentialEquations solvers](@id benchmark-ODE-solvers)

In this page, we will benchmark several solvers provided by [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/) for solving time evolution in hierarchical equations of motion approach.

```@example benchmark_ODE_solvers
using OrdinaryDiffEq ## or "using DifferentialEquations"
using BenchmarkTools
using HierarchicalEOM
HierarchicalEOM.versioninfo()
```
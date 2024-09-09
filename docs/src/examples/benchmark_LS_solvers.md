# [LinearSolve solvers](@id benchmark-LS-solvers)

In this page, we will benchmark several solvers provided by [LinearSolve.jl](https://docs.sciml.ai/LinearSolve/stable/) for solving SteadyState and spectrum in hierarchical equations of motion approach.

```@example benchmark_LS_solvers
using LinearSolve
using BenchmarkTools
using HierarchicalEOM
HierarchicalEOM.versioninfo()
```

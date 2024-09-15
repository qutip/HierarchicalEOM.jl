# # [DifferentialEquations solvers](@id benchmark-ODE-solvers)
# 
# In this page, we will benchmark several solvers provided by [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/) for solving time evolution in hierarchical equations of motion approach.  

using OrdinaryDiffEq ## or "using DifferentialEquations" 
using BenchmarkTools
using HierarchicalEOM
HierarchicalEOM.versioninfo()

# Here, we use the example of [driven systems and dynamical decoupling](@ref exp-dynamical-decoupling):

Γ = 0.0005
W = 0.005
kT = 0.05
N = 3
tier = 6
amp = 0.50
delay = 20
tlist = 0:0.4:400

σz = sigmaz()
σx = sigmax()
H0 = 0.0 * σz
ψ0 = (basis(2, 0) + basis(2, 1)) / √2
params = (V = amp, Δ = delay, σx = σx)

function pulse(V, Δ, t)
    τ = 0.5 * π / V
    period = τ + Δ

    if (t % period) < τ
        return V
    else
        return 0
    end
end

function H_D(t, p::NamedTuple)
    return pulse(p.V, p.Δ, t) * p.σx
end

bath = Boson_DrudeLorentz_Pade(σz, Γ, W, kT, N)
M = M_Boson(H0, tier, bath);

# ## ODE Solver List
# (click [here](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/) to see the full solver list provided by `DifferentialEquations.jl`)
#   
# For any extra solver options, we can add it in the function [`HEOMsolve`](@ref) with keyword arguments. These keyword arguments will be directly pass to the solvers in `DifferentialEquations`
# (click [here](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/) to see the documentation for the common solver options)
#   
# Furthermore, since we are using [`BenchmarkTools`](https://juliaci.github.io/BenchmarkTools.jl/stable/) (`@benchmark`) in the following, we set `verbose = false` to disable the output message.
#   
# ### DP5 (Default solver)
# Dormand-Prince's 5/4 Runge-Kutta method. (free 4th order interpolant)
@benchmark HEOMsolve(M, ψ0, tlist; H_t = H_D, params = params, abstol = 1e-12, reltol = 1e-12, verbose = false)

# ### RK4
# The canonical Runge-Kutta Order 4 method. Uses a defect control for adaptive stepping using maximum error over the whole interval.
@benchmark HEOMsolve(M, ψ0, tlist; H_t = H_D, params = params, solver = RK4(), abstol = 1e-12, reltol = 1e-12, verbose = false)

# ### Tsit5
# Tsitouras 5/4 Runge-Kutta method. (free 4th order interpolant).
@benchmark HEOMsolve(M, ψ0, tlist; H_t = H_D, params = params, solver = Tsit5(), abstol = 1e-12, reltol = 1e-12, verbose = false)

# ### Vern7
# Verner's “Most Efficient” 7/6 Runge-Kutta method. (lazy 7th order interpolant).
@benchmark HEOMsolve(M, ψ0, tlist; H_t = H_D, params = params, solver = Vern7(), abstol = 1e-12, reltol = 1e-12, verbose = false)

# ### Vern9
# Verner's “Most Efficient” 9/8 Runge-Kutta method. (lazy 9th order interpolant)
@benchmark HEOMsolve(M, ψ0, tlist; H_t = H_D, params = params, solver = Vern9(), abstol = 1e-12, reltol = 1e-12, verbose = false)

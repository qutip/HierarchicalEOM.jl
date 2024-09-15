# # Quick Start
# 
# ### Content
#  - [Import HierarchicalEOM.jl](#Import-HierarchicalEOM.jl)
#  - [System and Bath](#System-and-Bath)
#  - [HEOM Liouvillian superoperator](#HEOM-Liouvillian-superoperator)
#  - [Time Evolution](#Time-Evolution)
#  - [Stationary State](#Stationary-State)
#  - [Reduced Density Operator](#Reduced-Density-Operator)
#  - [Expectation Value](#Expectation-Value)
#  - [Multiple Baths](#Multiple-Baths)
# ### Import HierarchicalEOM.jl
# Here are the functions in `HierarchicalEOM.jl` that we will use in this tutorial (Quick Start):

import HierarchicalEOM
import HierarchicalEOM: Boson_DrudeLorentz_Pade, M_Boson, HEOMsolve, getRho, BosonBath

# Note that you can also type `using HierarchicalEOM` to import everything you need in `HierarchicalEOM.jl`.
# To check the versions of dependencies of `HierarchicalEOM.jl`, run the following function

HierarchicalEOM.versioninfo()

# ### System and Bath
# Let us consider a simple two-level system coupled to a Drude-Lorentz bosonic bath. The system Hamiltonian, ``H_{sys}``, and the bath spectral density, ``J_D``, are

# ```math
# H_{sys}=\frac{\epsilon \sigma_z}{2} + \frac{\Delta \sigma_x}{2} ~~\text{and}
# ```  

# ```math
# J_{D}(\omega)=\frac{2\lambda W\omega}{W^2+\omega^2},
# ```
# #### System Hamiltonian and initial state
# You must construct system hamiltonian, initial state, and coupling operators by [`QuantumToolbox`](https://github.com/qutip/QuantumToolbox.jl) framework. It provides many useful functions to create arbitrary quantum states and operators which can be combined in all the expected ways.

import QuantumToolbox: Qobj, sigmaz, sigmax, basis, ket2dm, expect, steadystate

## The system Hamiltonian
ϵ = 0.5 # energy of 2-level system
Δ = 1.0 # tunneling term

Hsys = 0.5 * ϵ * sigmaz() + 0.5 * Δ * sigmax()

## System initial state
ρ0 = ket2dm(basis(2, 0));

## Define the operators that measure the populations of the two system states:
P00 = ket2dm(basis(2, 0))
P11 = ket2dm(basis(2, 1))

## Define the operator that measures the 0, 1 element of density matrix
## (corresponding to coherence):
P01 = basis(2, 0) * basis(2, 1)'

# #### Bath Properties
# Now, we demonstrate how to describe the bath using the built-in implementation of ``J_D(\omega)`` under Pade expansion by calling [`Boson_DrudeLorentz_Pade`](@ref)

λ = 0.1  # coupling strength
W = 0.5  # band-width (cut-off frequency)
kT = 0.5  # the product of the Boltzmann constant k and the absolute temperature T

Q = sigmaz() # system-bath coupling operator

N = 2 # Number of expansion terms to retain:

## Padé expansion:
bath = Boson_DrudeLorentz_Pade(Q, λ, W, kT, N)

# For other different expansions of the different spectral density correlation functions, please refer to [Bosonic Bath](@ref doc-Bosonic-Bath) and [Fermionic Bath](@ref doc-Fermionic-Bath).

# ### HEOM Liouvillian superoperator
# For bosonic bath, we can construct the HEOM Liouvillian superoperator matrix by calling [`M_Boson`](@ref)

tier = 5 # maximum tier of hierarchy
L = M_Boson(Hsys, tier, bath)

# To learn more about the HEOM Liouvillian superoperator matrix (including other types: `M_Fermion`, `M_Boson_Fermion`), please refer to [HEOMLS Matrices](@ref doc-HEOMLS-Matrix).

# ### Time Evolution
# Next, we can calculate the time evolution for the entire auxiliary density operators (ADOs) by calling [`HEOMsolve`](@ref)
tlist = 0:0.2:50
sol = HEOMsolve(L, ρ0, tlist; e_ops = [P00, P11, P01])

# To learn more about `HEOMsolve`, please refer to [Time Evolution](@ref doc-Time-Evolution).

# ### Stationary State
# We can also solve the stationary state of the auxiliary density operators (ADOs) by calling [`steadystate`](@ref).
ados_steady = steadystate(L)

# To learn more about `steadystate`, please refer to [Stationary State](@ref doc-Stationary-State).

# ### Reduced Density Operator
# To obtain the reduced density operator, one can either access the first element of auxiliary density operator (`ADOs`) or call [`getRho`](@ref):

## reduce density operator in the final time (`end`) of the evolution
ados_list = sol.ados
ρ = ados_list[end][1]  # index `1` represents the reduced density operator
ρ = getRho(ados_list[end])

## reduce density operator in stationary state
ρ = ados_steady[1]
ρ = getRho(ados_steady)

# One of the great features of `HierarchicalEOM.jl` is that we allow users to not only considering the density operator of the reduced
# state but also easily take high-order terms into account without struggling in finding the indices (see [Auxiliary Density Operators](@ref doc-ADOs) and [Hierarchy Dictionary](@ref doc-Hierarchy-Dictionary) for more details).

# ### Expectation Value
# We can now compare the results obtained from `HEOMsolve` and `steadystate`:
## for time evolution
p00_e = real(sol.expect[1,:]) # P00 is the 1st element in e_ops
p01_e = real(sol.expect[3,:]); # P01 is the 3rd element in e_ops

## for steady state
p00_s = expect(P00, ados_steady)
p01_s = expect(P01, ados_steady);

# ### Plot the results
using Plots, LaTeXStrings

plot(tlist, p00_e, label = L"\textrm{P}_{00}", linecolor = :blue, linestyle = :solid, linewidth = 3, grid = false)
plot!(tlist, p01_e, label = L"\textrm{P}_{01}", linecolor = :red, linestyle = :solid, linewidth = 3)
plot!(
    tlist,
    ones(length(tlist)) .* p00_s,
    label = L"\textrm{P}_{00} \textrm{(Steady State)}",
    linecolor = :blue,
    linestyle = :dash,
    linewidth = 3,
)
plot!(
    tlist,
    ones(length(tlist)) .* p01_s,
    label = L"\textrm{P}_{01} \textrm{(Steady State)}",
    linecolor = :red,
    linestyle = :dash,
    linewidth = 3,
)

xlabel!("time")
ylabel!("Population")

# ### Multiple Baths
# `HierarchicalEOM.jl` also supports for system to interact with multiple baths.  

# All you need to do is to provide a list of baths instead of a single bath

## The system Hamiltonian
Hsys = Qobj([
    0.25 1.50 2.50
    1.50 0.75 3.50
    2.50 3.50 1.25
])

## System initial state
ρ0 = ket2dm(basis(3, 0));

## Projector for each system state:
P00 = ket2dm(basis(3, 0))
P11 = ket2dm(basis(3, 1))
P22 = ket2dm(basis(3, 2));

## Construct one bath for each system state:
## note that `BosonBath[]` make the list created in type: Vector{BosonBath}
baths = BosonBath[]
for i in 0:2
    ## system-bath coupling operator: |i><i|
    Q = ket2dm(basis(3, i))
    push!(baths, Boson_DrudeLorentz_Pade(Q, λ, W, kT, N))
end

L = M_Boson(Hsys, tier, baths)

tlist = 0:0.025:5
sol = HEOMsolve(L, ρ0, tlist; e_ops = [P00, P11, P22])

## calculate population for each system state:
p0 = real(sol.expect[1,:])
p1 = real(sol.expect[2,:])
p2 = real(sol.expect[3,:])

plot(tlist, p0, linewidth = 3, linecolor = "blue", label = L"P_0", grid = false)
plot!(tlist, p1, linewidth = 3, linecolor = "orange", label = L"P_1")
plot!(tlist, p2, linewidth = 3, linecolor = :green, label = L"P_2")
xlabel!("time")
ylabel!("Population")

# Note that this example can also be found in [qutip documentation](https://qutip.org/docs/latest/guide/heom/bosonic.html).

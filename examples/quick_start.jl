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
import HierarchicalEOM: Boson_DrudeLorentz_Pade, M_Boson, evolution, SteadyState, getRho, BosonBath, Expect

# Note that you can also type `using HierarchicalEOM` to import everything you need in `HierarchicalEOM.jl`.
# To check the versions of dependencies of `HierarchicalEOM.jl` , run the following function

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
# You can construct system hamiltonian, initial state, and coupling operators by standard julia bulit-in types: `Vector`, `SparseVector`, `Matrix`, `SparseMatrix`.

# Moreover, it is also convenient to use [`QuantumOptics`](https://qojulia.org/) or [`QuantumToolbox`](https://github.com/albertomercurio/QuantumToolbox.jl) framework. They both provide many useful functions to create arbitrary quantum states and operators which can be combined in all the expected ways.

# Note that `HierarchicalEOM.jl` only accept standard julia bulit-in types. If you use `QuantumOptics` or `QuantumToolbox` to construct the operators, remember that the matrix (or vector) in standard type are stored in `.data` field of the objects, i.e., `op.data`.  
# Thus, you should take it as the inputs of `HierarchicalEOM.jl`.

# !!! compat "Extension for QuantumOptics.jl"
#     `HierarchicalEOM.jl` provides an extension to support `QuantumOptics`-type object, but this feature requires `Julia 1.9+` and `HierarchicalEOM 0.3+`. See [here](@ref doc-ext-QuantumOptics) for more details.

# !!! compat "Extension for QuantumToolbox.jl"
#     `HierarchicalEOM.jl` provides an extension to support `QuantumToolbox`-type object, but this feature requires `Julia 1.9+` and `HierarchicalEOM 1.4+`. See [here](@ref doc-ext-QuantumToolbox) for more details.

# We demonstrate this tutorial by `QuantumOptics`:

import QuantumOptics: SpinBasis, sigmaz, sigmax, ⊗, Ket, Bra, dm

basis = SpinBasis(1 // 2)

## The system Hamiltonian
ϵ = 0.5 # energy of 2-level system
Δ = 1.0 # tunneling term

Hsys = 0.5 * ϵ * sigmaz(basis) + 0.5 * Δ * sigmax(basis)

## System initial state
ρ0 = dm(Ket(basis, [1, 0]));

# #### Bath Properties
# Now, we demonstrate how to describe the bath using the built-in implementation of ``J_D(\omega)`` under Pade expansion by calling [`Boson_DrudeLorentz_Pade`](@ref)

λ = 0.1  # coupling strength
W = 0.5  # band-width (cut-off frequency)
kT = 0.5  # the product of the Boltzmann constant k and the absolute temperature T

Q = sigmaz(basis) # system-bath coupling operator

N = 2 # Number of expansion terms to retain:

## Padé expansion:
## Remember to give the operator in Standard matrix (AbstractMatrix) type
## That is, if using QuantumOptics package, give system coupling operator as Q.data
bath = Boson_DrudeLorentz_Pade(Q.data, λ, W, kT, N)

# For other different expansions of the different spectral density correlation functions, please refer to [Bosonic Bath](@ref doc-Bosonic-Bath) and [Fermionic Bath](@ref doc-Fermionic-Bath).

# ### HEOM Liouvillian superoperator
# For bosonic bath, we can construct the HEOM Liouvillian superoperator matrix by calling [`M_Boson`](@ref)

## maximum tier of hierarchy
tier = 5

## Remember to give the operator in Standard matrix (AbstractMatrix) type
## That is, if using QuantumOptics package, give system hamiltonian as Hsys.data
L = M_Boson(Hsys.data, tier, bath)

# To learn more about the HEOM Liouvillian superoperator matrix (including other types: `M_Fermion`, `M_Boson_Fermion`), please refer to [HEOMLS Matrices](@ref doc-HEOMLS-Matrix).

# ### Time Evolution
# Next, we can calculate the time evolution for the entire auxiliary density operators (ADOs) by calling [`evolution`](@ref)
tlist = 0:0.2:50
ados_list = evolution(L, ρ0.data, tlist);

# To learn more about `evolution`, please refer to [Time Evolution](@ref doc-Time-Evolution).

# ### Stationary State
# We can also solve the stationary state of the auxiliary density operators (ADOs) by calling [`SteadyState`](@ref).
ados_steady = SteadyState(L)

# To learn more about `SteadyState`, please refer to [Stationary State](@ref doc-Stationary-State).

# ### Reduced Density Operator
# To obtain the reduced density operator, one can either access the first element of auxiliary density operator (`ADOs`) or call [`getRho`](@ref):

## reduce density operator in third time step of the evolution
ρ = ados_list[3][1]
ρ = getRho(ados_list[3])

## reduce density operator in stationary state
ρ = ados_steady[1]
ρ = getRho(ados_steady);

# One of the great features of `HierarchicalEOM.jl` is that we allow users to not only considering the density operator of the reduced
# state but also easily take high-order terms into account without struggling in finding the indices (see [Auxiliary Density Operators](@ref doc-ADOs) and [Hierarchy Dictionary](@ref doc-Hierarchy-Dictionary) for more details).

# ### Expectation Value
# We can now compare the results obtained from `evolution` and `SteadyState`:

## Define the operators that measure the populations of the two
## system states:
P00 = Ket(basis, [1, 0]) ⊗ Bra(basis, [1, 0])
P11 = Ket(basis, [0, 1]) ⊗ Bra(basis, [0, 1])

## Define the operator that measures the 0, 1 element of density matrix
## (corresponding to coherence):
P01 = Ket(basis, [1, 0]) ⊗ Bra(basis, [0, 1])

## for steady state
p00_s = Expect(P00.data, ados_steady)
p01_s = Expect(P01.data, ados_steady)

## for time evolution
p00_e = Expect(P00.data, ados_list)
p01_e = Expect(P01.data, ados_list);

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

# Note that, for the following, we use the built-in linear algebra in Julia (instead of `QuantumOptics.jl`) to construct the operators

## The system Hamiltonian
Hsys = [
    0.25 1.50 2.50
    1.50 0.75 3.50
    2.50 3.50 1.25
]

## System initial state
ρ0 = [
    1 0 0
    0 0 0
    0 0 0
];

## Construct one bath for each system state:
## note that `BosonBath[]` make the list created in type: Vector{BosonBath}
baths = BosonBath[]
for i in 1:3
    ## system-bath coupling operator: |i><i|
    Q = zeros(3, 3)
    Q[i, i] = 1

    push!(baths, Boson_DrudeLorentz_Pade(Q, λ, W, kT, N))
end

L = M_Boson(Hsys, tier, baths)

tlist = 0:0.025:5
ados_list = evolution(L, ρ0, tlist)

## Projector for each system state:
P00 = [1 0 0; 0 0 0; 0 0 0]
P11 = [0 0 0; 0 1 0; 0 0 0]
P22 = [0 0 0; 0 0 0; 0 0 1]

## calculate population for each system state:
p0 = Expect(P00, ados_list)
p1 = Expect(P11, ados_list)
p2 = Expect(P22, ados_list)

plot(tlist, p0, linewidth = 3, linecolor = "blue", label = L"P_0", grid = false)
plot!(tlist, p1, linewidth = 3, linecolor = "orange", label = L"P_1")
plot!(tlist, p2, linewidth = 3, linecolor = :green, label = L"P_2")
xlabel!("time")
ylabel!("Population")

# Note that this example can also be found in [qutip documentation](https://qutip.org/docs/latest/guide/heom/bosonic.html).

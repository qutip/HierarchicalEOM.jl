# [Extension for QuantumOptics.jl](@id doc-ext-QuantumOptics)

This is an extension to support [`QuantumOptics`](https://qojulia.org/)-type of operators (matrices)

!!! compat "Compat"
    The described feature requires `Julia 1.9+`.

When users construct the operators (such as system Hamiltonian, coupling operators, initial states, etc) by [`QuantumOptics.jl`](https://qojulia.org/) and take those matrices as inputs of `HierarchicalEOM.jl`, it will extract the `data` (matrix) part in `QuantumOptics.AbstractOperator`. Therefore, the type of `data` in `AbstractHEOMLSMatrix` still remain as `SparseMatrixCSC{ComplexF64, Int64}` without the basis information stored in `QuantumOptics.AbstractOperator`.

Although it doesn't store the basis information from `QuantumOptics.AbstractOperator`, it still supports calculating the expectation value ([`Expect`](@ref)) when the observable is given in the type of `QuantumOptics.AbstractOperator` (as long as the size of the matrix is equal).

Furthermore, it provides an extra method to re-construct reduced density operator with the type of `QuantumOptics.AbstractOperator`. Basically, it copies the basis information from another given `QuantumOptics.AbstractOperator` (as shown in the example below). With this functionality, one can again use the other functions provided in [`QuantumOptics.jl`](https://qojulia.org/).

The extension will be automatically loaded if user imports the package `QuantumOptics.jl` :
```@example
using QuantumOptics
using HierarchicalEOM

basis = SpinBasis(1//2)
I2 = identityoperator(basis)
d1 = sigmam(basis) ⊗      I2
d2 =      I2       ⊗ sigmam(basis)

# system Hamiltonian
Hsys = -3 * d1' * d1 - 5 * d2' * d2

# system initial state
ρ0 = dm(Ket(basis, [0, 1])) ⊗ dm(Ket(basis, [0, 1]))

# construct bath
bath1 = Fermion_Lorentz_Pade(d1, 0.01, 0, 10, 0.05, 2)
bath2 = Fermion_Lorentz_Pade(d2, 0.01, 0, 10, 0.05, 2)
baths = [bath1, bath2]

# construct HEOMLS matrix
tier = 2
L = M_Fermion(Hsys, tier, baths)
print(L)

# solving time evolution
tlist = 0:1:10
ados_list = evolution(L, ρ0, tlist)

# expectation value
exp_val = Expect(d1' * d1, ados_list)

# re-construct reduced density operator 
# with the type of QuantumOptics.AbstractOperator
ρ_normal  = ados_list[end][1] # reduced density operator ([1]) at t = 10 ([end]) 
ρ_QO_type = Operator(ρ_normal, Hsys) # basis information is copied from Hsys

# can use the functions provided in QuantumOptics.jl such as partial trace:
ptrace(ρ_QO_type, [1])
```
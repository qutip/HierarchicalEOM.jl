# [Extension for QuPhys.jl](@id doc-ext-QuPhys)

This is an extension to support [`QuPhys`](https://github.com/albertomercurio/QuPhys.jl)-type of operators (matrices)

!!! compat "Compat"
    The described feature requires `Julia 1.9+`.

When users construct the operators (such as system Hamiltonian, coupling operators, initial states, etc) by [`QuPhys.jl`](https://github.com/albertomercurio/QuPhys.jl) and take those matrices as inputs of `HierarchicalEOM.jl`, it will extract the `data` (matrix) part in `QuPhys.QuantumObject`. Therefore, the type of `data` in `AbstractHEOMLSMatrix` still remain as `SparseMatrixCSC{ComplexF64, Int64}` without the `type` and `dims` information stored in `QuPhys.QuantumObject`.

Although it doesn't store extra information (`type` and `dims`) from `QuPhys.QuantumObject`, it still supports calculating the expectation value ([`Expect`](@ref)) when the observable is given in the type of `QuPhys.QuantumObject` (as long as the size of the matrix is equal).

Furthermore, it provides an extra method to re-construct reduced density operator with the type of `QuPhys.QuantumObject`. Basically, it copies the `type` and `dims` information from another given `QuPhys.QuantumObject` (as shown in the example below). With this functionality, one can again use the other functions provided in [`QuPhys.jl`](https://github.com/albertomercurio/QuPhys.jl).

The extension will be automatically loaded if user imports the package `QuPhys.jl` :
```julia
using QuPhys
using HierarchicalEOM

d1 = sigmam() ⊗   eye(2)
d2 =   eye(2) ⊗ sigmam()

# system Hamiltonian
Hsys = -3 * d1' * d1 - 5 * d2' * d2

# system initial state
ρ0 = ket2dm(tensor(basis(2, 0), basis(2, 0)))

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
ρ_normal  = ados_list[end][1] # reduced density operator ([1]) at t = 10 ([end]) 

# with the type of QuPhys.QuantumObject (either the following two methods)
ρ_QO_type = Qobj(ρ_normal, Hsys)          # basis information is copied from Hsys
ρ_QO_type = QuantumObject(ρ_normal, Hsys) # basis information is copied from Hsys

# can use the functions provided in QuantumOptics.jl such as partial trace:
ptrace(ρ_QO_type, [1])
```
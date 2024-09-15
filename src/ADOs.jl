@doc raw"""
    struct ADOs
The Auxiliary Density Operators for HEOM model.

# Fields
- `data` : the vectorized auxiliary density operators
- `dims` : the dimension list of the coupling operator (should be equal to the system dims).
- `N` : the number of auxiliary density operators
- `parity`: the parity label (`EVEN` or `ODD`).

# Methods
One can obtain the density matrix for specific index (`idx`) by calling : `ados[idx]`.
`HierarchicalEOM.jl` also supports the following calls (methods) :
```julia
length(ados);  # returns the total number of `ADOs`
ados[1:idx];   # returns a vector which contains the `ADO` (in matrix form) from index `1` to `idx`
ados[1:end];   # returns a vector which contains all the `ADO` (in matrix form)
ados[:];       # returns a vector which contains all the `ADO` (in matrix form)
for rho in ados  # iteration
    # do something
end
```
"""
struct ADOs
    data::SparseVector{ComplexF64,Int64}
    dims::SVector
    N::Int
    parity::AbstractParity
end

# these functions are for forward compatibility
ADOs(data::SparseVector{ComplexF64,Int64}, dim::Int, N::Int, parity::AbstractParity) = ADOs(data, [dim], N, parity)
ADOs(data::SparseVector{ComplexF64,Int64}, dims::AbstractVector, N::Int, parity::AbstractParity) = ADOs(data, SVector{length(dims),Int}(dims), N, parity)

@doc raw"""
    ADOs(V, N, parity)
Gernerate the object of auxiliary density operators for HEOM model.

# Parameters
- `V::AbstractVector` : the vectorized auxiliary density operators
- `N::Int` : the number of auxiliary density operators.
- `parity::AbstractParity` : the parity label (`EVEN` or `ODD`). Default to `EVEN`.
"""
function ADOs(V::AbstractVector, N::Int, parity::AbstractParity = EVEN)
    # check the dimension of V
    d = size(V, 1)
    dim = √(d / N)
    if isinteger(dim)
        return ADOs(sparsevec(V), SVector{1,Int}(Int(dim)), N, parity)
    else
        error("The dimension of vector is not consistent with the ADOs number \"N\".")
    end
end

@doc raw"""
    ADOs(ρ, N, parity)
Gernerate the object of auxiliary density operators for HEOM model.

# Parameters
- `ρ` : the reduced density operator
- `N::Int` : the number of auxiliary density operators.
- `parity::AbstractParity` : the parity label (`EVEN` or `ODD`). Default to `EVEN`.
"""
function ADOs(ρ::QuantumObject, N::Int = 1, parity::AbstractParity = EVEN)
    _ρ = sparsevec(ket2dm(ρ).data)
    return ADOs(sparsevec(_ρ.nzind, _ρ.nzval, N * length(_ρ)), ρ.dims, N, parity)
end
ADOs(ρ, N::Int = 1, parity::AbstractParity = EVEN) =
    error("HierarchicalEOM doesn't support input `ρ` with type : $(typeof(ρ))")

checkbounds(A::ADOs, i::Int) =
    ((i > A.N) || (i < 1)) ? error("Attempt to access $(A.N)-element ADOs at index [$(i)]") : nothing

@doc raw"""
    length(A::ADOs)
Returns the total number of the Auxiliary Density Operators (ADOs)
"""
length(A::ADOs) = A.N

@doc raw"""
    eltype(A::ADOs)
Returns the elements' type of the Auxiliary Density Operators (ADOs)
"""
eltype(A::ADOs) = eltype(A.data)

lastindex(A::ADOs) = length(A)

function getindex(A::ADOs, i::Int)
    checkbounds(A, i)

    D = prod(A.dims)
    sup_dim = D^2
    back = sup_dim * i
    return QuantumObject(reshape(A.data[(back-sup_dim+1):back], D, D), Operator, A.dims)
end

function getindex(A::ADOs, r::UnitRange{Int})
    checkbounds(A, r[1])
    checkbounds(A, r[end])

    result = []
    D = prod(A.dims)
    sup_dim = D^2
    for i in r
        back = sup_dim * i
        push!(result, QuantumObject(reshape(A.data[(back-sup_dim+1):back], D, D), Operator, A.dims))
    end
    return result
end
getindex(A::ADOs, ::Colon) = getindex(A, 1:lastindex(A))

iterate(A::ADOs, state::Int = 1) = state > length(A) ? nothing : (A[state], state + 1)

show(io::IO, A::ADOs) = print(io, "$(A.N) Auxiliary Density Operators with $(A.parity) and (system) dims = $(A.dims)\n")
show(io::IO, m::MIME"text/plain", A::ADOs) = show(io, A)

@doc raw"""
    getRho(ados)
Return the density matrix of the reduced state (system) from a given auxiliary density operators

# Parameters
- `ados::ADOs` : the auxiliary density operators for HEOM model

# Returns
- `ρ::QuantumObject` : The density matrix of the reduced state
"""
function getRho(ados::ADOs)
    D = prod(ados.dims)
    return QuantumObject(reshape(ados.data[1:(D^2)], D, D), Operator, ados.dims)
end

@doc raw"""
    getADO(ados, idx)
Return the auxiliary density operator with a specific index from auxiliary density operators

This function equals to calling : `ados[idx]`.

# Parameters
- `ados::ADOs` : the auxiliary density operators for HEOM model
- `idx::Int` : the index of the auxiliary density operator

# Returns
- `ρ_idx::QuantumObject` : The auxiliary density operator
"""
getADO(ados::ADOs, idx::Int) = ados[idx]

@doc raw"""
    expect(op, ados; take_real=true)
Return the expectation value of the operator `op` for the reduced density operator in the given `ados`, namely
```math
\textrm{Tr}\left[ O \rho \right],
```
where ``O`` is the operator and ``\rho`` is the reduced density operator in the given ADOs.

# Parameters
- `op` : the operator ``O`` to take the expectation value
- `ados::ADOs` : the auxiliary density operators for HEOM model
- `take_real::Bool` : whether to automatically take the real part of the trace or not. Default to `true`

# Returns
- `exp_val` : The expectation value
"""
function expect(op, ados::ADOs; take_real::Bool = true)
    if typeof(op) <: HEOMSuperOp
        _check_sys_dim_and_ADOs_num(op, ados)
        exp_val = dot(transpose(_Tr(ados.dims, ados.N)), (SparseMatrixCSC(op) * ados).data)
    else
        _op = HandleMatrixType(op, ados.dims, "op (observable)"; type = Operator)
        exp_val = tr(_op.data * getRho(ados).data)
    end

    if take_real
        return real(exp_val)
    else
        return exp_val
    end
end

@doc raw"""
    expect(op, ados_list; take_real=true)
Return a list of expectation values of the operator `op` corresponds to the reduced density operators in the given `ados_list`, namely
```math
\textrm{Tr}\left[ O \rho \right],
```
where ``O`` is the operator and ``\rho`` is the reduced density operator in one of the `ADOs` from `ados_list`.

# Parameters
- `op` : the operator ``O`` to take the expectation value
- `ados_list::Vector{ADOs}` : the list of auxiliary density operators for HEOM model
- `take_real::Bool` : whether to automatically take the real part of the trace or not. Default to `true`

# Returns
- `exp_val` : The expectation value
"""
function expect(op, ados_list::Vector{ADOs}; take_real::Bool = true)
    dims = ados_list[1].dims
    N = ados_list[1].N
    for i in 2:length(ados_list)
        _check_sys_dim_and_ADOs_num(ados_list[1], ados_list[i])
    end

    if typeof(op) <: HEOMSuperOp
        _check_sys_dim_and_ADOs_num(op, ados_list[1])
        _op = op
    else
        _op = HEOMSuperOp(op, EVEN, dims, N, "L")
    end
    tr_op = transpose(_Tr(dims, N)) * SparseMatrixCSC(_op).data

    exp_val = [dot(tr_op, ados.data) for ados in ados_list]

    if take_real
        return real.(exp_val)
    else
        return exp_val
    end
end

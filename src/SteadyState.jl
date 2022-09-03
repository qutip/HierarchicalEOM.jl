"""
# `steadystate(M::AbstractHEOMMatrix; [solver, SOLVEROptions...])`
Solve the steady state of the given Heom matrix.

## Parameters
- `M::AbstractHEOMMatrix` : the matrix given from HEOM model (the parity should be either \":none (boson)\" or \":even (fermion)\".)
- `solver` : solver in package `LinearSolve.jl`. Default to `KLUFactorization()`.
- `SOLVEROptions` : extra options for solver 

## Returns
- `ADOs` : The auxiliary density operators of the steady state.
"""
function steadystate(M::AbstractHEOMMatrix; solver=KLUFactorization(), SOLVEROptions...)
    # check parity
    if (M.parity != :even) || (M.parity != :none)
        error("The parity of M should be either \":none (bonson)\" or \":even (fermion)\".")
    end    

    A = copy(M.data)
    S, = size(A)
    A[1,1:S] .= 0
    
    # sparse(row_idx, col_idx, values, row_dims, col_dims)
    A += sparse(fill(1, M.dim), [(n - 1) * (M.dim + 1) + 1 for n in 1:(M.dim)], fill(1, M.dim), S, S)
    
    b = sparsevec([1], [1. + 0.0im], s)
    
    # solving x where A * x = b
    sol = solve(LinearProblem(A, Vector(b)), solver, SOLVEROptions...)
    
    return ADOs(sol.u, M.dim, M.Nb, M.Nf)
end
"""
    Steadystate(M; solver, verbose, SOLVEROptions...)
Solve the steady state of the auxiliary density operators.

# Parameters
- `M::AbstractHEOMMatrix` : the matrix given from HEOM model, where the parity should be either `:none` (boson) or `:even` (fermion).
- `solver` : solver in package `LinearSolve.jl`. Default to `UMFPACKFactorization()`.
- `verbose::Bool` : To display verbose output and progress bar during the process or not. Defaults to `true`.
- `SOLVEROptions` : extra options for solver 

For more details about solvers and extra options, please refer to [`LinearSolve.jl`](http://linearsolve.sciml.ai/stable/)

# Returns
- `::ADOs` : The steady state of auxiliary density operators.
"""
function Steadystate(M::AbstractHEOMMatrix; solver=UMFPACKFactorization(), verbose::Bool=true, SOLVEROptions...)
    # check parity
    if (M.parity != :even) && (M.parity != :none)
        error("The parity of M should be either \":none\" (bonson) or \":even\" (fermion).")
    end    

    A = copy(M.data)
    S, = size(A)
    A[1,1:S] .= 0
    
    # sparse(row_idx, col_idx, values, row_dims, col_dims)
    A += sparse(fill(1, M.dim), [(n - 1) * (M.dim + 1) + 1 for n in 1:(M.dim)], fill(1, M.dim), S, S)
    
    b = sparsevec([1], [1. + 0.0im], S)
    
    # solving x where A * x = b
    if verbose
        print("Solving steady state for auxiliary density operators...")
        flush(stdout)
    end
    sol = solve(LinearProblem(A, Vector(b)), solver, SOLVEROptions...)
    if verbose
        println("[DONE]")
        flush(stdout)
    end
    
    return ADOs(sol.u, M.dim, M.Nb, M.Nf)
end
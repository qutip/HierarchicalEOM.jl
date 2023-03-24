@doc raw"""
    SteadyState(M; solver, verbose, SOLVEROptions...)
Solve the steady state of the auxiliary density operators based on `LinearSolve.jl` (i.e., solving ``x`` where ``A \times x = b``).

# Parameters
- `M::AbstractHEOMMatrix` : the matrix given from HEOM model, where the parity should be either `:none` (boson) or `:even` (fermion).
- `solver` : solver in package `LinearSolve.jl`. Default to `UMFPACKFactorization()`.
- `verbose::Bool` : To display verbose output and progress bar during the process or not. Defaults to `true`.
- `SOLVEROptions` : extra options for solver 

For more details about solvers and extra options, please refer to [`LinearSolve.jl`](http://linearsolve.sciml.ai/stable/)

# Returns
- `::ADOs` : The steady state of auxiliary density operators.
"""
@noinline function SteadyState(M::AbstractHEOMMatrix; solver=UMFPACKFactorization(), verbose::Bool=true, SOLVEROptions...)
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
    
    return ADOs(sol.u, M.dim, M.N)
end

# func. for solving ODE
function _hierarchy!(dρ, ρ, L, t)
    @inbounds dρ .= L * ρ
end

@doc raw"""
    SteadyState(M, ρ0; solver, reltol, abstol, maxiters, save_everystep, verbose, SOLVEROptions...)
Solve the steady state of the auxiliary density operators based on time evolution (ordinary differential equations)
with initial state is given in the type of density-matrix (`ρ0`).

# Parameters
- `M::AbstractHEOMMatrix` : the matrix given from HEOM model
- `ρ0` : system initial state (density matrix)
- `solver` : The ODE solvers in package `DifferentialEquations.jl`. Default to `FBDF(autodiff=false)`.
- `reltol::Real` : Relative tolerance in adaptive timestepping. Default to `1.0e-6`.
- `abstol::Real` : Absolute tolerance in adaptive timestepping. Default to `1.0e-8`.
- `maxiters::Real` : Maximum number of iterations before stopping. Default to `1e5`.
- `save_everystep::Bool` : Saves the result at every step. Defaults to `false`.
- `verbose::Bool` : To display verbose output and progress bar during the process or not. Defaults to `true`.
- `SOLVEROptions` : extra options for solver

For more details about solvers and extra options, please refer to [`DifferentialEquations.jl`](https://diffeq.sciml.ai/stable/)

# Returns
- `::ADOs` : The steady state of auxiliary density operators.
"""
function SteadyState(
        M::AbstractHEOMMatrix, 
        ρ0;
        solver = FBDF(autodiff=false),
        reltol::Real = 1.0e-6,
        abstol::Real = 1.0e-8,
        maxiters::Real = 1e5,
        save_everystep::Bool=false,
        verbose::Bool = true,
        SOLVEROptions...
    )

    if !isValidMatrixType(ρ0, M.dim)
        error("Invalid matrix \"ρ0\".")
    end

    # vectorize initial state
    ρ1   = sparse(sparsevec(ρ0))
    ados = ADOs(sparsevec(ρ1.nzind, ρ1.nzval, M.N * M.sup_dim), M.N)
    
    return SteadyState(M, ados;
        solver = solver,
        reltol = reltol,
        abstol = abstol,
        maxiters = maxiters,
        save_everystep = save_everystep,
        verbose = verbose,
        SOLVEROptions...
    )
end

@doc raw"""
    SteadyState(M, ados; solver, reltol, abstol, maxiters, save_everystep, verbose, SOLVEROptions...)
Solve the steady state of the auxiliary density operators based on time evolution (ordinary differential equations)
with initial state is given in the type of `ADOs`.

# Parameters
- `M::AbstractHEOMMatrix` : the matrix given from HEOM model
- `ados::ADOs` : initial auxiliary density operators
- `solver` : The ODE solvers in package `DifferentialEquations.jl`. Default to `FBDF(autodiff=false)`.
- `reltol::Real` : Relative tolerance in adaptive timestepping. Default to `1.0e-6`.
- `abstol::Real` : Absolute tolerance in adaptive timestepping. Default to `1.0e-8`.
- `maxiters::Real` : Maximum number of iterations before stopping. Default to `1e5`.
- `save_everystep::Bool` : Saves the result at every step. Defaults to `false`.
- `verbose::Bool` : To display verbose output and progress bar during the process or not. Defaults to `true`.
- `SOLVEROptions` : extra options for solver

For more details about solvers and extra options, please refer to [`DifferentialEquations.jl`](https://diffeq.sciml.ai/stable/)

# Returns
- `::ADOs` : The steady state of auxiliary density operators.
"""
@noinline function SteadyState(
        M::AbstractHEOMMatrix, 
        ados::ADOs;
        solver = FBDF(autodiff=false),
        reltol = 1.0e-6,
        abstol = 1.0e-8,
        maxiters = 1e5,
        save_everystep::Bool = false,
        verbose::Bool = true,
        SOLVEROptions...
    )
    
    # check parity
    if M.parity == :odd
        error("The parity of M should be either \":none\" (bonson) or \":even\" (fermion).")
    end

    if (M.dim != ados.dim)
        error("The system dimension between M and ados are not consistent.")
    end

    if (M.N != ados.N)
        error("The number N between M and ados are not consistent.")
    end

    # setup ode function
    jac_p = undef
    try
        if solver.linsolve == nothing
            jac_p = SparseMatrixCSC{ComplexF64, Int64}
        else
            S, = size(M)
            jac_p = spzeros(ComplexF64, S, S)
        end
    catch
        jac_p = SparseMatrixCSC{ComplexF64, Int64}
    end
    hierarchy = ODEFunction(_hierarchy!; jac_prototype = jac_p)

    # solving steady state of the ODE problem
    if verbose
        print("Solving steady state for auxiliary density operators...")
        flush(stdout)
    end
    sol = solve(
        SteadyStateProblem(hierarchy, Vector(ados.data), M.data), 
        DynamicSS(solver; abstol = abstol, reltol = reltol);
        maxiters = maxiters,
        save_everystep = save_everystep,
        SOLVEROptions...
    )

    if verbose
        println("[DONE]")
        flush(stdout)
    end

    return ADOs(sol.u, M.dim, M.N)
end
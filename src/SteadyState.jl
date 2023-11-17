@doc raw"""
    SteadyState(M; solver, verbose, SOLVEROptions...)
Solve the steady state of the auxiliary density operators based on `LinearSolve.jl` (i.e., solving ``x`` where ``A \times x = b``).

# Parameters
- `M::AbstractHEOMLSMatrix` : the matrix given from HEOM model, where the parity should be `EVEN`.
- `solver` : solver in package `LinearSolve.jl`. Default to `UMFPACKFactorization()`.
- `verbose::Bool` : To display verbose output and progress bar during the process or not. Defaults to `true`.
- `SOLVEROptions` : extra options for solver 

For more details about solvers and extra options, please refer to [`LinearSolve.jl`](http://linearsolve.sciml.ai/stable/)

# Returns
- `::ADOs` : The steady state of auxiliary density operators.
"""
@noinline function SteadyState(M::AbstractHEOMLSMatrix; solver=UMFPACKFactorization(), verbose::Bool=true, SOLVEROptions...)
    # check parity
    if typeof(M.parity) != EvenParity
        error("The parity of M should be \"EVEN\".")
    end    

    S = size(M, 1)
    A = _HandleSteadyStateMatrix(typeof(M.data), M, S)
    b = sparsevec([1], [1. + 0.0im], S)
    
    # solving x where A * x = b
    if verbose
        print("Solving steady state for auxiliary density operators...")
        flush(stdout)
    end
    cache = init(LinearProblem(A, _HandleVectorType(typeof(M.data), b)), solver, SOLVEROptions...)
    sol = solve!(cache)
    if verbose
        println("[DONE]")
        flush(stdout)
    end
    
    return ADOs(_HandleVectorType(sol.u, false), M.dim, M.N, M.parity)
end

@doc raw"""
    SteadyState(M, ρ0; solver, reltol, abstol, maxiters, save_everystep, verbose, SOLVEROptions...)
Solve the steady state of the auxiliary density operators based on time evolution (ordinary differential equations)
with initial state is given in the type of density-matrix (`ρ0`).

# Parameters
- `M::AbstractHEOMLSMatrix` : the matrix given from HEOM model, where the parity should be `EVEN`.
- `ρ0` : system initial state (density matrix)
- `solver` : The ODE solvers in package `DifferentialEquations.jl`. Default to `DP5()`.
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
        M::AbstractHEOMLSMatrix, 
        ρ0;
        solver = DP5(),
        reltol::Real = 1.0e-6,
        abstol::Real = 1.0e-8,
        maxiters::Real = 1e5,
        save_everystep::Bool=false,
        verbose::Bool = true,
        SOLVEROptions...
    )
    return SteadyState(
        M, 
        ADOs(ρ0, M.N, M.parity);
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
- `M::AbstractHEOMLSMatrix` : the matrix given from HEOM model, where the parity should be `EVEN`.
- `ados::ADOs` : initial auxiliary density operators
- `solver` : The ODE solvers in package `DifferentialEquations.jl`. Default to `DP5()`.
- `reltol::Real` : Relative tolerance in adaptive timestepping. Default to `1.0e-3`.
- `abstol::Real` : Absolute tolerance in adaptive timestepping. Default to `1.0e-6`.
- `maxiters::Real` : Maximum number of iterations before stopping. Default to `1e5`.
- `save_everystep::Bool` : Saves the result at every step. Defaults to `false`.
- `verbose::Bool` : To display verbose output and progress bar during the process or not. Defaults to `true`.
- `SOLVEROptions` : extra options for solver

For more details about solvers and extra options, please refer to [`DifferentialEquations.jl`](https://diffeq.sciml.ai/stable/)

# Returns
- `::ADOs` : The steady state of auxiliary density operators.
"""
@noinline function SteadyState(
        M::AbstractHEOMLSMatrix, 
        ados::ADOs;
        solver = DP5(),
        reltol = 1.0e-6,
        abstol = 1.0e-8,
        maxiters = 1e5,
        save_everystep::Bool = false,
        verbose::Bool = true,
        SOLVEROptions...
    )
    
    _check_sys_dim_and_ADOs_num(M, ados)
    _check_parity(M, ados)

    # problem: dρ(t)/dt = L * ρ(t)
    L = MatrixOperator(M.data)
    ElType = eltype(M)
    tspan  = Tuple{ElType, ElType}((0, Inf))
    prob   = ODEProblem(L, _HandleVectorType(typeof(M.data), ados.data), tspan)

    # solving steady state of the ODE problem
    if verbose
        print("Solving steady state for auxiliary density operators...")
        flush(stdout)
    end
    sol = solve(
        SteadyStateProblem(prob), 
        DynamicSS(solver; abstol = _HandleFloatType(eltype(M), abstol), reltol = _HandleFloatType(eltype(M), reltol));
        maxiters = maxiters,
        save_everystep = save_everystep,
        SOLVEROptions...
    )

    if verbose
        println("[DONE]")
        flush(stdout)
    end

    return ADOs(_HandleVectorType(sol.u, false), M.dim, M.N, M.parity)
end
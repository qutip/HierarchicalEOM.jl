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
        print("Solving steady state for ADOs by linear-solve method...")
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
    SteadyState(M, ρ0, tspan; solver, termination_condition, verbose, SOLVEROptions...)
Solve the steady state of the auxiliary density operators based on time evolution (ordinary differential equations; `SteadyStateDiffEq.jl`) with initial state is given in the type of density-matrix (`ρ0`).

# Parameters
- `M::AbstractHEOMLSMatrix` : the matrix given from HEOM model, where the parity should be `EVEN`.
- `ρ0` : system initial state (density matrix)
- `tspan::Number` : the time limit to find stationary state. Default to `Inf`
- `solver` : The ODE solvers in package `DifferentialEquations.jl`. Default to `DP5()`.
- `termination_condition` : The stationary state terminate condition in `DiffEqBase.jl`. Default to `NormTerminationMode()`.
- `verbose::Bool` : To display verbose output and progress bar during the process or not. Defaults to `true`.
- `SOLVEROptions` : extra options for solver

For more details about solvers, termination condition, and extra options, please refer to [`DifferentialEquations.jl`](https://diffeq.sciml.ai/stable/)

# Returns
- `::ADOs` : The steady state of auxiliary density operators.
"""
function SteadyState(
        M::AbstractHEOMLSMatrix, 
        ρ0,
        tspan::Number = Inf;
        solver = DP5(),
        termination_condition = NormTerminationMode(),
        verbose::Bool = true,
        SOLVEROptions...
    )
    return SteadyState(
        M, 
        ADOs(ρ0, M.N, M.parity),
        tspan;
        solver = solver,
        termination_condition = termination_condition,
        verbose = verbose,
        SOLVEROptions...
    )
end

@doc raw"""
    SteadyState(M, ados, tspan; solver, termination_condition, verbose, SOLVEROptions...)
Solve the steady state of the auxiliary density operators based on time evolution (ordinary differential equations; `SteadyStateDiffEq.jl`) with initial state is given in the type of `ADOs`.

# Parameters
- `M::AbstractHEOMLSMatrix` : the matrix given from HEOM model, where the parity should be `EVEN`.
- `ados::ADOs` : initial auxiliary density operators
- `tspan::Number` : the time limit to find stationary state. Default to `Inf`
- `solver` : The ODE solvers in package `DifferentialEquations.jl`. Default to `DP5()`.
- `termination_condition` : The stationary state terminate condition in `DiffEqBase.jl`. Default to `NormTerminationMode()`.
- `verbose::Bool` : To display verbose output and progress bar during the process or not. Defaults to `true`.
- `SOLVEROptions` : extra options for solver

For more details about solvers, termination condition, and extra options, please refer to [`DifferentialEquations.jl`](https://diffeq.sciml.ai/stable/)

# Returns
- `::ADOs` : The steady state of auxiliary density operators.
"""
@noinline function SteadyState(
        M::AbstractHEOMLSMatrix, 
        ados::ADOs,
        tspan::Number = Inf;
        solver = DP5(),
        termination_condition = NormTerminationMode(),
        verbose::Bool = true,
        SOLVEROptions...
    )
    
    _check_sys_dim_and_ADOs_num(M, ados)
    _check_parity(M, ados)

    # problem: dρ(t)/dt = L * ρ(t)
    L = MatrixOperator(M.data)
    prob = SteadyStateProblem(L, _HandleVectorType(typeof(M.data), ados.data))

    # solving steady state of the ODE problem
    if verbose
        print("Solving steady state for ADOs by Ordinary Differential Equations method...")
        flush(stdout)
    end
    sol = solve(
        prob, 
        DynamicSS(solver, _HandleFloatType(eltype(M), tspan));
        termination_condition = termination_condition,
        SOLVEROptions...
    )

    if verbose
        println("[DONE]")
        flush(stdout)
    end

    return ADOs(_HandleVectorType(sol.u, false), M.dim, M.N, M.parity)
end
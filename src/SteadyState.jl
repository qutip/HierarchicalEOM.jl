@doc raw"""
    SteadyState(M; solver, verbose, SOLVEROptions...)
Solve the steady state of the auxiliary density operators based on `LinearSolve.jl` (i.e., solving ``x`` where ``A \times x = b``).

# Parameters
- `M::AbstractHEOMLSMatrix` : the matrix given from HEOM model, where the parity should be `EVEN`.
- `solver::SciMLLinearSolveAlgorithm` : solver in package `LinearSolve.jl`. Default to `UMFPACKFactorization()`.
- `verbose::Bool` : To display verbose output and progress bar during the process or not. Defaults to `true`.
- `SOLVEROptions` : extra options for solver 

# Notes
- For more details about `solver` and `SOLVEROptions`, please refer to [`LinearSolve.jl`](http://linearsolve.sciml.ai/stable/)

# Returns
- `::ADOs` : The steady state of auxiliary density operators.
"""
@noinline function SteadyState(
    M::AbstractHEOMLSMatrix;
    solver::SciMLLinearSolveAlgorithm = UMFPACKFactorization(),
    verbose::Bool = true,
    SOLVEROptions...,
)
    # check parity
    if typeof(M.parity) != EvenParity
        error("The parity of M should be \"EVEN\".")
    end

    S = size(M, 1)
    A = _HandleSteadyStateMatrix(M, S)
    b = sparsevec([1], [1.0 + 0.0im], S)

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

    return ADOs(_HandleVectorType(sol.u, false), M.dims, M.N, M.parity)
end

@doc raw"""
    SteadyState(M, ρ0, tspan; solver, verbose, SOLVEROptions...)
Solve the steady state of the auxiliary density operators based on time evolution (`OrdinaryDiffEq.jl`) with initial state is given in the type of density-matrix (`ρ0`).

# Parameters
- `M::AbstractHEOMLSMatrix` : the matrix given from HEOM model, where the parity should be `EVEN`.
- `ρ0::QuantumObject` : system initial state (density matrix)
- `tspan::Number` : the time limit to find stationary state. Default to `Inf`
- `solver::OrdinaryDiffEqAlgorithm` : The ODE solvers in package `DifferentialEquations.jl`. Default to `DP5()`.
- `reltol::Real` : Relative tolerance in adaptive timestepping. Default to `1.0e-8`.
- `abstol::Real` : Absolute tolerance in adaptive timestepping. Default to `1.0e-10`.
- `save_everystep::Bool` : Saves the result at every step. Defaults to `false`.
- `verbose::Bool` : To display verbose output and progress bar during the process or not. Defaults to `true`.
- `SOLVEROptions` : extra options for solver

# Notes
- For more details about `solver` please refer to [`DifferentialEquations.jl` (ODE Solvers)](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/)
- For more details about `SOLVEROptions` please refer to [`DifferentialEquations.jl` (Keyword Arguments)](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/)

# Returns
- `::ADOs` : The steady state of auxiliary density operators.
"""
function SteadyState(
    M::AbstractHEOMLSMatrix,
    ρ0::QuantumObject,
    tspan::Number = Inf;
    solver::OrdinaryDiffEqAlgorithm = DP5(),
    reltol::Real = 1.0e-8,
    abstol::Real = 1.0e-10,
    save_everystep::Bool = false,
    verbose::Bool = true,
    SOLVEROptions...,
)
    return SteadyState(
        M,
        ADOs(ρ0, M.N, M.parity),
        tspan;
        solver = solver,
        reltol = reltol,
        abstol = abstol,
        save_everystep = save_everystep,
        verbose = verbose,
        SOLVEROptions...,
    )
end

@doc raw"""
    SteadyState(M, ados, tspan; solver, verbose, SOLVEROptions...)
Solve the steady state of the auxiliary density operators based on time evolution (`OrdinaryDiffEq.jl`) with initial state is given in the type of `ADOs`.

# Parameters
- `M::AbstractHEOMLSMatrix` : the matrix given from HEOM model, where the parity should be `EVEN`.
- `ados::ADOs` : initial auxiliary density operators
- `tspan::Number` : the time limit to find stationary state. Default to `Inf`
- `solver::OrdinaryDiffEqAlgorithm` : The ODE solvers in package `DifferentialEquations.jl`. Default to `DP5()`.
- `reltol::Real` : Relative tolerance in adaptive timestepping. Default to `1.0e-8`.
- `abstol::Real` : Absolute tolerance in adaptive timestepping. Default to `1.0e-10`.
- `save_everystep::Bool` : Saves the result at every step. Defaults to `false`.
- `verbose::Bool` : To display verbose output and progress bar during the process or not. Defaults to `true`.
- `SOLVEROptions` : extra options for solver

# Notes
- For more details about `solver` please refer to [`DifferentialEquations.jl` (ODE Solvers)](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/)
- For more details about `SOLVEROptions` please refer to [`DifferentialEquations.jl` (Keyword Arguments)](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/)

# Returns
- `::ADOs` : The steady state of auxiliary density operators.
"""
@noinline function SteadyState(
    M::AbstractHEOMLSMatrix,
    ados::ADOs,
    tspan::Number = Inf;
    solver::OrdinaryDiffEqAlgorithm = DP5(),
    reltol::Real = 1.0e-8,
    abstol::Real = 1.0e-10,
    save_everystep::Bool = false,
    verbose::Bool = true,
    SOLVEROptions...,
)
    _check_sys_dim_and_ADOs_num(M, ados)
    _check_parity(M, ados)

    ElType = eltype(M)
    Tspan = (_HandleFloatType(ElType, 0), _HandleFloatType(ElType, tspan))
    RelTol = _HandleFloatType(ElType, reltol)
    AbsTol = _HandleFloatType(ElType, abstol)

    # problem: dρ(t)/dt = L * ρ(t)
    prob = ODEProblem{true}(MatrixOperator(M.data), _HandleVectorType(typeof(M.data), ados.data), Tspan)

    # solving steady state of the ODE problem
    if verbose
        println("Solving steady state for ADOs by Ordinary Differential Equations method...")
        flush(stdout)
    end
    sol = solve(
        prob,
        solver;
        callback = TerminateSteadyState(AbsTol, RelTol, _ss_condition),
        reltol = RelTol,
        abstol = AbsTol,
        save_everystep = save_everystep,
        SOLVEROptions...,
    )

    if verbose
        println("Last timepoint t = $(sol.t[end])\n[DONE]")
        flush(stdout)
    end

    return ADOs(_HandleVectorType(sol.u[end], false), M.dims, M.N, M.parity)
end

function _ss_condition(integrator, abstol, reltol, min_t)
    # this condition is same as DiffEqBase.NormTerminationMode

    du_dt = (integrator.u - integrator.uprev) / integrator.dt
    norm_du_dt = norm(du_dt)
    if (norm_du_dt <= reltol * norm(du_dt + integrator.u)) || (norm_du_dt <= abstol)
        return true
    else
        return false
    end
end

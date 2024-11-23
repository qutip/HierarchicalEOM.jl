@doc raw"""
    steadystate(M::AbstractHEOMLSMatrix; solver, verbose, SOLVEROptions...)
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
function QuantumToolbox.steadystate(
    M::AbstractHEOMLSMatrix{<:MatrixOperator};
    solver::SciMLLinearSolveAlgorithm = UMFPACKFactorization(),
    verbose::Bool = true,
    SOLVEROptions...,
)
    # check parity
    if typeof(M.parity) != EvenParity
        error("The parity of M should be \"EVEN\".")
    end

    A = _HandleSteadyStateMatrix(M)
    b = sparsevec([1], [1.0 + 0.0im], size(M, 1))

    # solving x where A * x = b
    if verbose
        print("Solving steady state for ADOs by linear-solve method...")
        flush(stdout)
    end
    cache = init(LinearProblem(A, _HandleVectorType(M, b)), solver, SOLVEROptions...)
    sol = solve!(cache)
    if verbose
        println("[DONE]")
        flush(stdout)
    end

    return ADOs(Vector{ComplexF64}(sol.u), M.dims, M.N, M.parity)
end

@doc raw"""
    steadystate(M::AbstractHEOMLSMatrix, ρ0, tspan; solver, verbose, SOLVEROptions...)
Solve the steady state of the auxiliary density operators based on time evolution (`OrdinaryDiffEq.jl`) with initial state is given in the type of density-matrix (`ρ0`).

# Parameters
- `M::AbstractHEOMLSMatrix` : the matrix given from HEOM model, where the parity should be `EVEN`.
- `ρ0::Union{QuantumObject,ADOs}` : system initial state (density matrix) or initial auxiliary density operators (`ADOs`)
- `tspan::Number` : the time limit to find stationary state. Default to `Inf`
- `solver::OrdinaryDiffEqAlgorithm` : The ODE solvers in package `DifferentialEquations.jl`. Default to `DP5()`.
- `verbose::Bool` : To display verbose output and progress bar during the process or not. Defaults to `true`.
- `SOLVEROptions` : extra options for solver

# Notes
- For more details about `solver` please refer to [`DifferentialEquations.jl` (ODE Solvers)](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/)
- For more details about `SOLVEROptions` please refer to [`DifferentialEquations.jl` (Keyword Arguments)](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/)

# Returns
- `::ADOs` : The steady state of auxiliary density operators.
"""
function QuantumToolbox.steadystate(
    M::AbstractHEOMLSMatrix{<:MatrixOperator},
    ρ0::T_state,
    tspan::Number = Inf;
    solver::OrdinaryDiffEqAlgorithm = DP5(),
    verbose::Bool = true,
    SOLVEROptions...,
) where {T_state<:Union{QuantumObject,ADOs}}

    # handle initial state
    ados = (T_state <: QuantumObject) ? ADOs(ρ0, M.N, M.parity) : ρ0
    _check_sys_dim_and_ADOs_num(M, ados)
    _check_parity(M, ados)
    u0 = _HandleVectorType(M, ados.data)

    Tspan = (0, tspan)

    kwargs = merge(DEFAULT_ODE_SOLVER_OPTIONS, SOLVEROptions)
    cb = TerminateSteadyState(kwargs.abstol, kwargs.reltol, _ss_condition)
    kwargs2 =
        haskey(kwargs, :callback) ? merge(kwargs, (callback = CallbackSet(kwargs.callback, cb),)) :
        merge(kwargs, (callback = cb,))

    # define ODE problem
    prob = ODEProblem{true,FullSpecialize}(M.data, u0, Tspan; kwargs2...)

    # solving steady state of the ODE problem
    if verbose
        println("Solving steady state for ADOs by Ordinary Differential Equations method...")
        flush(stdout)
    end
    sol = solve(prob, solver)

    if verbose
        println("Last timepoint t = $(sol.t[end])\n[DONE]")
        flush(stdout)
    end

    return ADOs(Vector{ComplexF64}(sol.u[end]), M.dims, M.N, M.parity)
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

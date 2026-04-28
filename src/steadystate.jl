@doc raw"""
    steadystate(M::AbstractHEOMLSMatrix; alg, verbose, kwargs...)
Solve the steady state of the auxiliary density operators based on `LinearSolve.jl` (i.e., solving ``x`` where ``A \times x = b``).

# Parameters
- `M::AbstractHEOMLSMatrix` : the matrix given from HEOM model, where the parity should be `EVEN`.
- `alg::SciMLLinearSolveAlgorithm` : The solving algorithm in package `LinearSolve.jl`. Default to `KrylovJL_GMRES(rtol = 1e-12, atol = 1e-14)`.
- `verbose::Bool` : To display verbose output or not. Defaults to `true`.
- `kwargs` : The keyword arguments for the `LinearProblem`

# Notes
- For more details about `alg`, `kwargs`, and `LinearProblem`, please refer to [`LinearSolve.jl`](http://linearsolve.sciml.ai/stable/).
- This function supports [lazy operators](@ref doc-Lazy-Operators) for memory-efficient calculations. When using lazy operators, the algorithm must be a matrix-free solver (e.g., Krylov-based methods like `KrylovJL_GMRES`) that does not require concretizing the matrix.

# Returns
- `::ADOs` : The steady state of auxiliary density operators.
"""
function QuantumToolbox.steadystate(
        M::AbstractHEOMLSMatrix;
        alg::SciMLLinearSolveAlgorithm = KrylovJL_GMRES(rtol = 1.0e-12, atol = 1.0e-14),
        verbose::Bool = true,
        kwargs...,
    )
    isconstant(M) || throw(ArgumentError("The HEOMLS matrix M must be time-independent to solve steadystate."))
    haskey(kwargs, :solver) &&
        error("The keyword argument `solver` for solving HEOM steadystate is deprecated, use `alg` instead.")

    # check parity
    if typeof(M.parity) != EvenParity
        error("The parity of M should be \"EVEN\".")
    end

    b = _HandleVectorType(M, sparsevec([1], [1.0 + 0.0im], size(M, 1)))
    A = _HandleSteadyStateMatrix(M, b)

    # solving x where A * x = b
    if verbose
        print("Solving steady state for ADOs by linear-solve method...")
        flush(stdout)
    end

    prob = LinearProblem{true}(A, b)
    sol = solve(prob, alg; kwargs...)
    if verbose
        println("[DONE]")
        flush(stdout)
    end

    return ADOs(Vector{ComplexF64}(sol.u), M.dimensions, M.N, M.parity)
end

@doc raw"""
    steadystate(M::AbstractHEOMLSMatrix, ρ0, tspan; alg, terminate_reltol, terminate_abstol, verbose, kwargs...)
Solve the steady state of the auxiliary density operators based on time evolution (`OrdinaryDiffEq.jl`) with initial state is given in the type of density-matrix (`ρ0`).

# Parameters
- `M::AbstractHEOMLSMatrix` : the matrix given from HEOM model, where the parity should be `EVEN`. Supports both full sparse matrices and lazy tensor product representations (constructed with `assemble = Val(:combine)`).
- `ρ0::Union{QuantumObject,ADOs}` : system initial state (density matrix) or initial auxiliary density operators (`ADOs`)
- `tspan::Number` : the time limit to find stationary state. Default to `Inf`
- `alg::AbstractODEAlgorithm` : The ODE algorithm in package `DifferentialEquations.jl`. Default to `DP5()`.
- `terminate_reltol` = The relative tolerance for stationary state terminate condition. Default to `1e-6`.
- `terminate_abstol` = The absolute tolerance for stationary state terminate condition. Default to `1e-8`.
- `verbose::Bool` : To display verbose output or not. Defaults to `true`.
- `kwargs` : The keyword arguments in `ODEProblem`

# Notes
- For more details about `alg` please refer to [`DifferentialEquations.jl` (ODE Solvers)](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/)
- For more details about `kwargs` please refer to [`DifferentialEquations.jl` (Keyword Arguments)](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/)
- This function supports [lazy operators](@ref doc-Lazy-Operators) for memory-efficient calculations.

# Returns
- `::ADOs` : The steady state of auxiliary density operators.
"""
function QuantumToolbox.steadystate(
        M::AbstractHEOMLSMatrix,
        ρ0::T_state,
        tspan::Number = Inf;
        alg::AbstractODEAlgorithm = DP5(),
        terminate_reltol::Real = 1e-6,
        terminate_abstol::Real = 1e-8,
        verbose::Bool = true,
        kwargs...,
    ) where {T_state <: Union{QuantumObject, ADOs}}
    isconstant(M) || throw(ArgumentError("The HEOMLS matrix M must be time-independent to solve steadystate."))
    haskey(kwargs, :solver) &&
        error("The keyword argument `solver` for solving HEOM steadystate is deprecated, use `alg` instead.")

    (typeof(M.parity) == EvenParity) || error("The parity of M should be \"EVEN\".")

    # handle initial state
    ados = (T_state <: QuantumObject) ? ADOs(ρ0, M.N, EVEN) : ρ0
    _check_sys_dim_and_ADOs_num(M, ados)
    _check_parity(M, ados)
    u0 = _HandleVectorType(M, ados.data)

    ftype = _float_type(M)
    Tspan = (ftype(0), ftype(tspan))
    ode_options = default_ode_solver_options(ftype)

    kwargs2 = merge(
        (
            abstol = ode_options.abstol,
            reltol = ode_options.reltol,
            save_everystep = false,
            saveat = ftype[],
        ),
        kwargs,
    )

    # handle callbacks
    _ss_condition = SteadyStateODECondition(similar(u0))
    cb = TerminateSteadyState(terminate_abstol, terminate_reltol, _ss_condition)
    kwargs3 =
        haskey(kwargs2, :callback) ? merge(kwargs2, (callback = CallbackSet(kwargs2.callback, cb),)) :
        merge(kwargs2, (callback = cb,))

    A = get_cached_HEOMLS_data(M.data, u0)

    # define ODE problem
    prob = ODEProblem{true, FullSpecialize}(A, u0, Tspan; kwargs3...)

    # solving steady state of the ODE problem
    if verbose
        println("Solving steady state for ADOs by Ordinary Differential Equations method...")
        flush(stdout)
    end
    sol = solve(prob, alg)
    if verbose
        println("Last time point t = $(sol.t[end])\n[DONE]")
        flush(stdout)
    end

    return ADOs(Vector{ComplexF64}(sol.u[end]), M.dimensions, M.N, M.parity)
end

export HEOMsolve, heomsolve, HEOMsolve_map, heomsolve_map, TimeEvolutionHEOMSol, HEOMsolveProblem

@doc raw"""
    struct TimeEvolutionHEOMSol

A structure storing the results and some information from solving time evolution of hierarchical equations of motion (HEOM).

# Fields (Attributes)
- `Btier` : The tier (cutoff level) for bosonic hierarchy
- `Ftier` : The tier (cutoff level) for fermionic hierarchy
- `times::AbstractVector`: The list of time points at which the expectation values are calculated during the evolution.
- `times_ados::AbstractVector`: The list of time points at which the [`ADOs`](@ref) are stored during the evolution.
- `ados::Vector{ADOs}`: The list of result [`ADOs`](@ref) corresponding to each time point in `times_ados`.
- `expect::Union{AbstractMatrix,Nothing}`: The expectation values corresponding to each time point in `times`.
- `retcode`: The return code from the solver.
- `alg`: The algorithm which is used during the solving process.
- `abstol::Real`: The absolute tolerance which is used during the solving process.
- `reltol::Real`: The relative tolerance which is used during the solving process.
"""
struct TimeEvolutionHEOMSol{
    TT1<:AbstractVector{<:Real},
    TT2<:AbstractVector{<:Real},
    TS<:Vector{ADOs},
    TE<:Union{Nothing,AbstractMatrix},
    RETT<:Union{Nothing,Enum},
    AlgT<:Union{Nothing,AbstractODEAlgorithm},
    TolT<:Union{Nothing,Real},
}
    Btier::Int
    Ftier::Int
    times::TT1
    times_ados::TT2
    ados::TS
    expect::TE
    retcode::RETT
    alg::AlgT
    abstol::TolT
    reltol::TolT
end

function Base.show(io::IO, sol::TimeEvolutionHEOMSol)
    print(io, "Solution of hierarchical EOM\n")
    print(io, "(return code: $(sol.retcode))\n")
    print(io, "----------------------------\n")
    print(io, "Btier = $(sol.Btier)\n")
    print(io, "Ftier = $(sol.Ftier)\n")
    print(io, "num_ados   = $(length(sol.ados))\n")
    if sol.expect isa Nothing
        print(io, "num_expect = 0\n")
    else
        print(io, "num_expect = $(size(sol.expect, 1))\n")
    end
    print(io, "ODE alg.: $(sol.alg)\n")
    print(io, "abstol = $(sol.abstol)\n")
    print(io, "reltol = $(sol.reltol)\n")
    return nothing
end

function _gen_ados_ode_vector(ρ::QuantumObject, M::AbstractHEOMLSMatrix)
    ados = ADOs(ρ, M.N, M.parity)
    return _HandleVectorType(M, ados.data)
end
function _gen_ados_ode_vector(ados::ADOs, M::AbstractHEOMLSMatrix)
    _check_sys_dim_and_ADOs_num(M, ados)
    _check_parity(M, ados)
    return _HandleVectorType(M, ados.data)
end

@doc raw"""
    HEOMsolveProblem(
        M::AbstractHEOMLSMatrix,
        ρ0::Union{QuantumObject,ADOs},
        tlist::AbstractVector;
        e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
        H_t::Union{Nothing,QuantumObjectEvolution} = nothing,
        params = NullParameters(),
        progress_bar::Union{Val,Bool} = Val(true),
        inplace::Union{Val,Bool} = Val(true),
        kwargs...,
    )

Generate the ODEProblem for the time evolution of auxiliary density operators.

# Parameters
- `M::AbstractHEOMLSMatrix` : the matrix given from HEOM model
- `ρ0::Union{QuantumObject,ADOs}` : system initial state (density matrix) or initial auxiliary density operators (`ADOs`)
- `tlist::AbstractVector` : Denote the specific time points to save the solution at, during the solving process.
- `e_ops::Union{Nothing,AbstractVector}`: List of operators for which to calculate expectation values.
- `alg::AbstractODEAlgorithm` : The solving algorithm in package `DifferentialEquations.jl`. Default to `DP5()`.
- `H_t::Union{Nothing,QuantumObjectEvolution}`: The time-dependent system Hamiltonian or Liouvillian. Default to `nothing`.
- `params`: Parameters to pass to the solver. This argument is usually expressed as a `NamedTuple` or `AbstractVector` of parameters. For more advanced usage, any custom struct can be used.
- `progress_bar::Union{Val,Bool}`: Whether to show the progress bar. Defaults to `Val(true)`. Using non-`Val` types might lead to type instabilities.
- `inplace::Union{Val,Bool}`: Whether to use the inplace version of the ODEProblem. Defaults to `Val(true)`.
- `kwargs` : The keyword arguments for the `ODEProblem`.

# Notes
- The [`ADOs`](@ref) will be saved depend on the keyword argument `saveat` in `kwargs`.
- If `e_ops` is specified, the default value of `saveat=[tlist[end]]` (only save the final `ADOs`), otherwise, `saveat=tlist` (saving the `ADOs` corresponding to `tlist`). You can also specify `e_ops` and `saveat` separately.
- The default tolerances in `kwargs` are given as `reltol=1e-6` and `abstol=1e-8`.
- For more details about `alg` please refer to [`DifferentialEquations.jl` (ODE Solvers)](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/)
- For more details about `kwargs` please refer to [`DifferentialEquations.jl` (Keyword Arguments)](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/)

# Returns

- `prob`: The `TimeEvolutionProblem` containing the `ODEProblem` for the time evolution of auxiliary density operators.
"""
function HEOMsolveProblem(
    M::AbstractHEOMLSMatrix,
    ρ0::T_state,
    tlist::AbstractVector;
    e_ops::Union{Nothing,AbstractVector,Tuple} = nothing,
    H_t::Union{Nothing,QuantumObjectEvolution} = nothing,
    params = NullParameters(),
    progress_bar::Union{Val,Bool} = Val(true),
    inplace::Union{Val,Bool} = Val(true),
    kwargs...,
) where {T_state<:Union{QuantumObject,ADOs}}
    haskey(kwargs, :save_idxs) &&
        throw(ArgumentError("The keyword argument \"save_idxs\" is not supported in HierarchicalEOM."))

    tlist = _check_tlist(tlist, _float_type(M))
    tspan = (tlist[1], tlist[end])

    # handle initial state
    u0 = _gen_ados_ode_vector(ρ0, M)

    # define ODE problem (L should be an AbstractSciMLOperator)
    L = cache_operator(_make_L(M, H_t), u0)
    kwargs2 = _merge_saveat(tlist, e_ops, DEFAULT_ODE_SOLVER_OPTIONS; kwargs...)
    kwargs3 = _merge_tstops(kwargs2, isconstant(L), tlist)
    kwargs4 = _generate_heom_kwargs(e_ops, makeVal(progress_bar), tlist, kwargs3, SaveFuncHEOMSolve, M)
    prob = ODEProblem{getVal(inplace),FullSpecialize}(L, u0, tspan, params; kwargs4...)

    return TimeEvolutionProblem(prob, tlist, ADOsType(), M.dimensions, (M = M,))
end

@doc raw"""
    HEOMsolve(M, ρ0, tlist; e_ops, alg, H_t, params, progress_bar, inplace, kwargs...)
    heomsolve(M, ρ0, tlist; e_ops, alg, H_t, params, progress_bar, inplace, kwargs...)

Solve the time evolution for auxiliary density operators based on ordinary differential equations.

# Parameters
- `M::AbstractHEOMLSMatrix` : the matrix given from HEOM model
- `ρ0::Union{QuantumObject,ADOs}` : system initial state (density matrix) or initial auxiliary density operators (`ADOs`)
- `tlist::AbstractVector` : Denote the specific time points to save the solution at, during the solving process.
- `e_ops::Union{Nothing,AbstractVector}`: List of operators for which to calculate expectation values.
- `alg::AbstractODEAlgorithm` : The solving algorithm in package `DifferentialEquations.jl`. Default to `DP5()`.
- `H_t::Union{Nothing,QuantumObjectEvolution}`: The time-dependent system Hamiltonian or Liouvillian. Default to `nothing`.
- `params`: Parameters to pass to the solver. This argument is usually expressed as a `NamedTuple` or `AbstractVector` of parameters. For more advanced usage, any custom struct can be used.
- `progress_bar::Union{Val,Bool}`: Whether to show the progress bar. Defaults to `Val(true)`. Using non-`Val` types might lead to type instabilities.
- `inplace::Union{Val,Bool}`: Whether to use the inplace version of the ODEProblem. Defaults to `Val(true)`.
- `kwargs` : The keyword arguments for the `ODEProblem`.

# Notes
- The [`ADOs`](@ref) will be saved depend on the keyword argument `saveat` in `kwargs`.
- If `e_ops` is specified, the default value of `saveat=[tlist[end]]` (only save the final `ADOs`), otherwise, `saveat=tlist` (saving the `ADOs` corresponding to `tlist`). You can also specify `e_ops` and `saveat` separately.
- The default tolerances in `kwargs` are given as `reltol=1e-6` and `abstol=1e-8`.
- For more details about `alg` please refer to [`DifferentialEquations.jl` (ODE Solvers)](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/)
- For more details about `kwargs` please refer to [`DifferentialEquations.jl` (Keyword Arguments)](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/)

# Returns
- sol::TimeEvolutionHEOMSol : The solution of the hierarchical EOM. See also [`TimeEvolutionHEOMSol`](@ref)

!!! note
    `heomsolve` is a synonym of `HEOMsolve`.
"""
function HEOMsolve(
    M::AbstractHEOMLSMatrix,
    ρ0::T_state,
    tlist::AbstractVector;
    e_ops::Union{Nothing,AbstractVector} = nothing,
    alg::AbstractODEAlgorithm = DP5(),
    H_t::Union{Nothing,QuantumObjectEvolution} = nothing,
    params = NullParameters(),
    progress_bar::Union{Val,Bool} = Val(true),
    inplace::Union{Val,Bool} = Val(true),
    kwargs...,
) where {T_state<:Union{QuantumObject,ADOs}}
    haskey(kwargs, :solver) && error("The keyword argument `solver` for HEOMsolve is deprecated, use `alg` instead.")
    haskey(kwargs, :verbose) &&
        error("The keyword argument `verbose` for HEOMsolve is deprecated, use `progress_bar` instead.")
    haskey(kwargs, :filename) && error("The keyword argument `filename` for HEOMsolve is deprecated.")

    # Move sensealg argument to solve for Enzyme.jl support.
    # TODO: Remove it when https://github.com/SciML/SciMLSensitivity.jl/issues/1225 is fixed.
    sensealg = get(kwargs, :sensealg, nothing)
    kwargs_filtered = isnothing(sensealg) ? kwargs : Base.structdiff((; kwargs...), (sensealg = sensealg,))

    prob = HEOMsolveProblem(
        M,
        ρ0,
        tlist;
        e_ops = e_ops,
        H_t = H_t,
        params = params,
        progress_bar = progress_bar,
        inplace = inplace,
        kwargs_filtered...,
    )

    # TODO: Remove it when https://github.com/SciML/SciMLSensitivity.jl/issues/1225 is fixed.
    if isnothing(sensealg)
        return HEOMsolve(prob, alg)
    else
        return HEOMsolve(prob, alg; sensealg = sensealg)
    end
end

function HEOMsolve(prob::TimeEvolutionProblem, alg::AbstractODEAlgorithm = DP5(); kwargs...)
    sol = solve(prob.prob, alg; kwargs...)

    return _gen_HEOMsolve_solution(sol, prob.times, prob.kwargs.M)
end

function _gen_HEOMsolve_solution(sol, times, M::AbstractHEOMLSMatrix)
    ADOs_list = map(ρvec -> ADOs(Vector{ComplexF64}(ρvec), M.dimensions, M.N, M.parity), sol.u)

    kwargs = NamedTuple(sol.prob.kwargs) # Convert to NamedTuple for Zygote.jl compatibility

    return TimeEvolutionHEOMSol(
        _getBtier(M),
        _getFtier(M),
        times,
        sol.t,
        ADOs_list,
        _get_expvals(sol, SaveFuncHEOMSolve),
        sol.retcode,
        sol.alg,
        kwargs.abstol,
        kwargs.reltol,
    )
end

const heomsolve = HEOMsolve # a synonym to align with qutip

function _generate_Eops(M::AbstractHEOMLSMatrix, e_ops)
    tr_e_ops = [
        # another adjoint will be applied in dot function in the HEOMsolveCallback
        _HandleTraceVectorType(M, adjoint(HEOMSuperOp(spre(op), EVEN, M.dimensions, M.N).data) * _Tr(M)) for op in e_ops
    ]
    return tr_e_ops
end

struct SaveFuncHEOMSolve{TE,PT<:Union{Nothing,Progress},IT,TEXPV<:Union{Nothing,AbstractMatrix}} <: AbstractSaveFunc
    tr_e_ops::TE
    progr::PT
    iter::IT
    expvals::TEXPV
end

(f::SaveFuncHEOMSolve)(u, t, integrator) = _save_func_heomsolve(u, integrator, f.tr_e_ops, f.progr, f.iter, f.expvals)
(f::SaveFuncHEOMSolve{Nothing})(u, t, integrator) = _save_func(integrator, f.progr) # Common for both mesolve and sesolve

function _save_func_heomsolve(u, integrator, tr_e_ops, progr, iter, expvals)
    _expect = op -> dot(op, u)
    @. expvals[:, iter[]] = _expect(tr_e_ops)
    iter[] += 1
    return _save_func(integrator, progr)
end

function _generate_heom_kwargs(
    e_ops,
    progress_bar,
    tlist,
    kwargs,
    method::Type{SaveFuncHEOMSolve},
    M::AbstractHEOMLSMatrix,
)
    tr_e_ops = e_ops isa Nothing ? nothing : _generate_Eops(M, e_ops)

    progr =
        getVal(progress_bar) ?
        Progress(
            length(tlist);
            enabled = getVal(progress_bar),
            desc = "[HEOMsolve (ODE)] ",
            QuantumToolbox.settings.ProgressMeterKWARGS...,
        ) : nothing

    expvals = e_ops isa Nothing ? nothing : Array{ComplexF64}(undef, length(e_ops), length(tlist))

    _save_func = method(tr_e_ops, progr, Ref(1), expvals)
    cb = FunctionCallingCallback(_save_func, funcat = tlist)

    return _merge_kwargs_with_callback(kwargs, cb)
end
_generate_heom_kwargs(
    e_ops::Nothing,
    progress_bar::Val{false},
    tlist,
    kwargs,
    method::Type{SaveFuncHEOMSolve},
    M::AbstractHEOMLSMatrix,
) = kwargs

_make_L(M::AbstractHEOMLSMatrix, H_t::Nothing) = M.data
function _make_L(M::AbstractHEOMLSMatrix, H_t::QuantumObjectEvolution)
    MType = _get_SciML_matrix_wrapper(M)
    L_t = HandleMatrixType(liouvillian(H_t), M.dimensions, "H_t"; type = SuperOperator())

    return M.data + _L_t_to_HEOMSuperOp(MType, L_t.data, M.N)
end

_L_t_to_HEOMSuperOp(MType::Type{<:AbstractSparseMatrix}, M::MatrixOperator, N::Int) =
    MatrixOperator(MType(kron(Eye(N), M.A)))
_L_t_to_HEOMSuperOp(MType::Type{<:AbstractSparseMatrix}, M::ScaledOperator, N::Int) =
    ScaledOperator(M.λ, _L_t_to_HEOMSuperOp(MType, M.L, N))
_L_t_to_HEOMSuperOp(MType::Type{<:AbstractSparseMatrix}, M::AddedOperator, N::Int) =
    AddedOperator(map(op -> _L_t_to_HEOMSuperOp(MType, op, N), M.ops))

@doc raw"""
    HEOMsolve_map(M, ρ0, tlist; e_ops, alg, ensemblealg, H_t, params, progress_bar, kwargs...)
    heomsolve_map(M, ρ0, tlist; e_ops, alg, ensemblealg, H_t, params, progress_bar, kwargs...)

Solve the time evolution for auxiliary density operators with multiple initial states and parameter sets using ensemble simulation.

This function computes the time evolution for all combinations (Cartesian product) of initial states and parameter sets, solving the HEOM (see [`HEOMsolve`](@ref)) for each combination in the ensemble.

# Arguments
- `M::AbstractHEOMLSMatrix` : the matrix given from HEOM model
- `ρ0::AbstractVector{<:Union{QuantumObject,ADOs}}` : system initial state(s) [the elements can be density matrix or `ADOs`].
- `tlist::AbstractVector` : Denote the specific time points to save the solution at, during the solving process.
- `e_ops::Union{Nothing,AbstractVector}`: List of operators for which to calculate expectation values.
- `alg::AbstractODEAlgorithm` : The solving algorithm in package `DifferentialEquations.jl`. Default to `DP5()`.
- `ensemblealg::EnsembleAlgorithm`: Ensemble algorithm to use for parallel computation. Default is `EnsembleThreads()`.
- `H_t::Union{Nothing,QuantumObjectEvolution}`: The time-dependent system Hamiltonian or Liouvillian. Default to `nothing`.
- `params::Union{NullParameters,Tuple}`: A `Tuple` of parameter sets. Each element should be an `AbstractVector` representing the sweep range for that parameter. The function will solve for all combinations of initial states and parameter sets.
- `progress_bar::Union{Val,Bool}`: Whether to show the progress bar. Default to `Val(true)`. Using non-`Val` types might lead to type instabilities.
- `kwargs` : The keyword arguments for the `ODEProblem`.

# Notes
- The function returns an array of solutions with dimensions matching the Cartesian product of initial states and parameter sets.
- If `ρ0` is a vector with length `m`, and `params = (p1, p2, ...)` where `p1` has length `n1`, `p2` has length `n2`, etc., the output will be of size `(m, n1, n2, ...)`.
- See [`HEOMsolve`](@ref) for more details.

# Returns
- An array of [`TimeEvolutionHEOMSol`](@ref) objects with dimensions `(length(ρ0), length(params[1]), length(params[2]), ...)`.

!!! note
    `heomsolve_map` is a synonym of `HEOMsolve_map`.
"""
function HEOMsolve_map(
    M::AbstractHEOMLSMatrix,
    ρ0::AbstractVector{<:T_state},
    tlist::AbstractVector;
    e_ops::Union{Nothing,AbstractVector} = nothing,
    alg::AbstractODEAlgorithm = DP5(),
    ensemblealg::EnsembleAlgorithm = EnsembleThreads(),
    H_t::Union{Nothing,QuantumObjectEvolution} = nothing,
    params::Union{NullParameters,Tuple} = NullParameters(),
    progress_bar::Union{Val,Bool} = Val(true),
    kwargs...,
) where {T_state<:Union{QuantumObject,ADOs}}

    # mapping initial states and parameters
    ados_iter = map(ρ -> _gen_ados_ode_vector(ρ, M), ρ0)
    if params isa NullParameters
        iter = collect(Iterators.product(ados_iter, [params])) |> vec # convert nx1 Matrix into Vector
    else
        iter = collect(Iterators.product(ados_iter, params...))
    end

    # we disable the progress bar of the HEOMsolveProblem because we use a global progress bar for all the trajectories
    prob = HEOMsolveProblem(
        M,
        first(ρ0),
        tlist;
        e_ops = e_ops,
        H_t = H_t,
        params = first(iter)[2:end],
        progress_bar = Val(false),
        kwargs...,
    )

    return HEOMsolve_map(prob, iter, alg, ensemblealg; progress_bar = progress_bar)
end
HEOMsolve_map(
    M::AbstractHEOMLSMatrix,
    ρ0::T_state,
    tlist::AbstractVector;
    kwargs...,
) where {T_state<:Union{QuantumObject,ADOs}} = HEOMsolve_map(M, [ρ0], tlist; kwargs...)

# this method is for advanced usage
# User can define their own iterator structure, prob_func and output_func
#   - `prob_func`: Function to use for generating the ODEProblem.
#   - `output_func`: a `Tuple` containing the `Function` to use for generating the output of a single trajectory, the (optional) `Progress` object, and the (optional) `RemoteChannel` object.
#
# Return: An array of TimeEvolutionSol objects with the size same as the given iter.
function HEOMsolve_map(
    prob::TimeEvolutionProblem{<:QuantumObjectType,<:AbstractDimensions,<:ODEProblem},
    iter::AbstractArray,
    alg::AbstractODEAlgorithm = DP5(),
    ensemblealg::EnsembleAlgorithm = EnsembleThreads();
    prob_func::Union{Function,Nothing} = nothing,
    output_func::Union{Tuple,Nothing} = nothing,
    progress_bar::Union{Val,Bool} = Val(true),
)
    # generate ensemble problem
    ntraj = length(iter)
    _prob_func = isnothing(prob_func) ? (prob, i, repeat) -> _se_me_map_prob_func(prob, i, repeat, iter) : prob_func
    _output_func =
        isnothing(output_func) ?
        _ensemble_dispatch_output_func(
            ensemblealg,
            progress_bar,
            ntraj,
            _standard_output_func;
            progr_desc = "[HEOMsolve_map] ",
        ) : output_func
    ens_prob = TimeEvolutionProblem(
        EnsembleProblem(prob.prob, prob_func = _prob_func, output_func = _output_func[1], safetycopy = false),
        prob.times,
        ADOsType(),
        prob.dimensions,
        (progr = _output_func[2], channel = _output_func[3]),
    )

    sol = _ensemble_dispatch_solve(ens_prob, alg, ensemblealg, ntraj)

    # handle solution and make it become an Array of TimeEvolutionHEOMSol
    sol_vec = [_gen_HEOMsolve_solution(sol[:, i], prob.times, prob.kwargs.M) for i in eachindex(sol)] # map is type unstable
    return reshape(sol_vec, size(iter))
end

const heomsolve_map = HEOMsolve_map # a synonym to align with QuantumToolbox.jl

module HierarchicalEOM

# Standard Julia libraries
using LinearAlgebra
using SparseArrays

import Base.Threads: @threads, nthreads, Channel
import Pkg

# Re-export QuantumToolbox
import Reexport: @reexport
@reexport using QuantumToolbox

# intrinsic QuantumToolbox functions
import QuantumToolboxCore:
    QuantumToolboxCore,
    _float_type,
    _complex_float_type,
    _spre,
    _spost,
    _sprepost,
    _liouvillian,
    _sum_lindblad_dissipators,
    _gen_dimensions,
    _get_dims_string,
    dimensions_to_dims,
    get_size,
    makeVal,
    getVal
import QuantumToolbox:
    _check_tlist,
    _save_func,
    _merge_saveat,
    _merge_tstops,
    _merge_kwargs_with_callback,
    _get_expvals,
    _se_me_map_prob_func,
    _standard_output_func,
    _ensemble_dispatch_output_func,
    _ensemble_dispatch_solve,
    TimeEvolutionProblem,
    AbstractSaveFunc,
    default_ode_solver_options,
    SteadyStateODECondition

# SciML packages (for OrdinaryDiffEq and LinearSolve)
import SciMLBase:
    SciMLBase,
    init,
    solve,
    solve!,
    u_modified!,
    ODEProblem,
    EnsembleProblem,
    EnsembleAlgorithm,
    EnsembleThreads,
    FullSpecialize,
    CallbackSet,
    NullParameters,
    AbstractODEAlgorithm
import SciMLOperators:
    SciMLOperators,
    AbstractSciMLOperator,
    DiagonalOperator,
    MatrixOperator,
    ScaledOperator,
    IdentityOperator,
    TensorProductOperator,
    AddedOperator,
    update_coefficients!,
    concretize
import OrdinaryDiffEqLowOrderRK: DP5
import DiffEqCallbacks: FunctionCallingCallback, TerminateSteadyState
import LinearSolve: LinearSolve, LinearProblem, needs_concrete_A, SciMLLinearSolveAlgorithm, KrylovJL_GMRES

# other dependencies (in alphabetical order)
import FastExpm: fastExpm
import FillArrays: Eye
import IncompleteLU: ilu
import ProgressMeter: Progress, next!

# Basic functions
include("HeomBase.jl")

# Bath
include("bath/BathBase.jl")
include("bath/BosonBath.jl")
include("bath/FermionBath.jl")
include("bath_correlation_functions/bath_correlation_func.jl")

# Parity and ADOs
include("Parity.jl")
include("ADOs.jl")

# HEOM Liouvillian superoperator (HEOMLS) matrices
include("heom_matrices/heom_matrix_base.jl")
include("heom_matrices/Nvec.jl")
include("heom_matrices/HierarchyDict.jl")
include("heom_matrices/M_S.jl")
include("heom_matrices/M_Boson.jl")
include("heom_matrices/M_Fermion.jl")
include("heom_matrices/M_Boson_Fermion.jl")

# Solvers
include("evolution.jl")
include("evolution_propagator.jl")
include("steadystate.jl")
include("power_spectrum.jl")
include("density_of_states.jl")
include("correlations.jl")

# deprecated functions
include("deprecated.jl")

function __init__()
    # register QuantumToolbox library and its dependencies
    if (HierarchicalEOM ∉ QuantumToolboxCore.QT_LIBRARIES)
        # use pushfirst! so that main API libraries are at the front of the registry (for better display order in versioninfo)
        pushfirst!(QuantumToolboxCore.QT_LIBRARIES, HierarchicalEOM)

        # dependencies
        m_list = Module[SciMLBase, SciMLOperators, LinearSolve]
        foreach(m_list) do m
            (m ∉ QuantumToolboxCore.DEP_PKGS) && push!(QuantumToolboxCore.DEP_PKGS, m)
        end
    end

    return nothing
end

end

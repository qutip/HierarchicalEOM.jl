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
import QuantumToolbox:
    _float_type,
    _complex_float_type,
    _check_tlist,
    _spre,
    _spost,
    _sprepost,
    _liouvillian,
    _sum_lindblad_dissipators,
    _gen_dimensions,
    _get_dims_string,
    dimensions_to_dims,
    _save_func,
    _merge_saveat,
    _merge_tstops,
    _merge_kwargs_with_callback,
    _get_expvals,
    _se_me_map_prob_func,
    _standard_output_func,
    _ensemble_dispatch_output_func,
    _ensemble_dispatch_solve,
    makeVal,
    getVal,
    TimeEvolutionProblem,
    AbstractSaveFunc,
    DEFAULT_ODE_SOLVER_OPTIONS,
    SteadyStateODECondition

# SciML packages (for OrdinaryDiffEq and LinearSolve)
import SciMLBase:
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
    MatrixOperator,
    ScaledOperator,
    AddedOperator,
    update_coefficients!,
    concretize
import OrdinaryDiffEqLowOrderRK: DP5
import DiffEqCallbacks: FunctionCallingCallback, TerminateSteadyState
import LinearSolve: LinearProblem, SciMLLinearSolveAlgorithm, KrylovJL_GMRES

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
include("steadystate.jl")
include("power_spectrum.jl")
include("density_of_states.jl")
include("correlations.jl")

# deprecated functions
include("deprecated.jl")

end

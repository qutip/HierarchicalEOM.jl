module HierarchicalEOM

# Re-export QuantumToolbox
import Reexport: @reexport
@reexport using QuantumToolbox

# intrinsic QuantumToolbox functions
import QuantumToolbox:
    _FType,
    _CType,
    _check_tlist,
    _spre,
    _spost,
    _sprepost,
    _liouvillian,
    _gen_dimensions,
    _get_dims_string,
    dimensions_to_dims,
    _merge_saveat,
    DEFAULT_ODE_SOLVER_OPTIONS

# SciML packages (for OrdinaryDiffEq and LinearSolve)
import SciMLBase: init, solve, solve!, u_modified!, ODEProblem, FullSpecialize, CallbackSet, NullParameters
import SciMLOperators:
    AbstractSciMLOperator, MatrixOperator, ScaledOperator, AddedOperator, update_coefficients!, concretize
import OrdinaryDiffEqCore: OrdinaryDiffEqAlgorithm
import OrdinaryDiffEqLowOrderRK: DP5
import DiffEqCallbacks: PresetTimeCallback, TerminateSteadyState
import LinearSolve: LinearProblem, SciMLLinearSolveAlgorithm, UMFPACKFactorization

# other dependencies (in alphabetical order)
import Base.Threads: @threads, threadid, nthreads, lock, unlock, SpinLock
import FastExpm: fastExpm
import JLD2: jldopen
import Pkg

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

# deprecated functions
include("deprecated.jl")

end

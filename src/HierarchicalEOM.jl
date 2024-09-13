module HierarchicalEOM

import Reexport: @reexport
@reexport using QuantumToolbox

# sub-module HeomBase for HierarchicalEOM
module HeomBase
    import Pkg
    import LinearAlgebra: BLAS, kron, I
    import SparseArrays: sparse, SparseVector, SparseMatrixCSC
    import StaticArraysCore: SVector
    import QuantumToolbox: QuantumObject, QuantumObjectType, Operator, SuperOperator

    export _Tr,
        PROGBAR_OPTIONS,
        AbstractHEOMLSMatrix,
        _check_sys_dim_and_ADOs_num,
        _check_parity,
        HandleMatrixType,
        _HandleFloatType,
        _HandleVectorType,
        _HandleSteadyStateMatrix,
        _HandleIdentityType

    include("HeomBase.jl")
end
import .HeomBase.versioninfo as versioninfo
import .HeomBase.print_logo as print_logo
@reexport import .HeomBase: AbstractHEOMLSMatrix

# sub-module Bath for HierarchicalEOM
module Bath
    using ..HeomBase
    import Base: show, length, getindex, lastindex, iterate, checkbounds
    import LinearAlgebra: ishermitian, eigvals, I
    import SparseArrays: SparseMatrixCSC
    import StaticArraysCore: SVector
    import QuantumToolbox: QuantumObject, _spre, _spost

    export AbstractBath,
        BosonBath,
        BosonBathRWA,
        FermionBath,
        Exponent,
        C,
        AbstractBosonBath,
        bosonReal,
        bosonImag,
        bosonRealImag,
        bosonAbsorb,
        bosonEmit,
        AbstractFermionBath,
        fermionAbsorb,
        fermionEmit,
        Boson_DrudeLorentz_Matsubara,
        Boson_DrudeLorentz_Pade,
        Fermion_Lorentz_Matsubara,
        Fermion_Lorentz_Pade

    include("Bath.jl")
    include("bath_correlation_functions/bath_correlation_func.jl")
end
@reexport using .Bath

# sub-module HeomAPI for HierarchicalEOM
module HeomAPI
    import Reexport: @reexport
    using ..HeomBase
    using ..Bath
    import Base:
        ==,
        !,
        +,
        -,
        *,
        show,
        length,
        size,
        getindex,
        keys,
        setindex!,
        lastindex,
        iterate,
        checkbounds,
        hash,
        copy,
        eltype
    import Base.Threads: @threads, threadid, nthreads, lock, unlock, SpinLock
    import LinearAlgebra: I, kron, tr, norm
    import SparseArrays: sparse, sparsevec, spzeros, SparseVector, SparseMatrixCSC
    import StaticArraysCore: SVector
    import QuantumToolbox:
        QuantumObject,
        Operator,
        SuperOperator,
        _spre,
        _spost,
        spre,
        spost,
        sprepost,
        expect,
        ket2dm,
        lindblad_dissipator
    import ProgressMeter: Progress, next!
    import FastExpm: fastExpm

    # solving time evolution and steady state
    import SciMLBase: init, solve, solve!, step!, ODEProblem
    import SciMLOperators: MatrixOperator
    import OrdinaryDiffEqCore: OrdinaryDiffEqAlgorithm
    import OrdinaryDiffEqLowOrderRK: DP5
    import DiffEqCallbacks: TerminateSteadyState
    import LinearSolve: LinearProblem, SciMLLinearSolveAlgorithm, UMFPACKFactorization
    import JLD2: jldopen

    export AbstractParity,
        OddParity,
        EvenParity,
        value,
        ODD,
        EVEN,
        ADOs,
        getRho,
        getADO,
        Nvec,
        AbstractHierarchyDict,
        HierarchyDict,
        MixHierarchyDict,
        getIndexEnsemble,
        HEOMSuperOp,
        M_S,
        M_Boson,
        M_Fermion,
        M_Boson_Fermion,
        Propagator,
        addBosonDissipator,
        addFermionDissipator,
        addTerminator,
        evolution,
        SteadyState,
        PowerSpectrum,
        DensityOfStates

    include("Parity.jl")
    include("ADOs.jl")

    include("heom_matrices/heom_matrix_base.jl")
    include("heom_matrices/Nvec.jl")
    include("heom_matrices/HierarchyDict.jl")
    include("heom_matrices/M_S.jl")
    include("heom_matrices/M_Boson.jl")
    include("heom_matrices/M_Fermion.jl")
    include("heom_matrices/M_Boson_Fermion.jl")

    include("evolution.jl")
    include("SteadyState.jl")
    include("power_spectrum.jl")
    include("density_of_states.jl")
end
@reexport using .HeomAPI

end

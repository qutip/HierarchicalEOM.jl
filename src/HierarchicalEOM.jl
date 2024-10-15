module HierarchicalEOM

import Reexport: @reexport
@reexport using QuantumToolbox

# sub-module HeomBase for HierarchicalEOM
module HeomBase
    import Pkg
    import LinearAlgebra: BLAS, kron, I
    import SparseArrays: sparse, SparseVector, SparseMatrixCSC, AbstractSparseMatrix
    import StaticArraysCore: SVector
    import QuantumToolbox: _FType, _CType, QuantumObject, QuantumObjectType, Operator, SuperOperator

    export _Tr,
        AbstractHEOMLSMatrix,
        _check_sys_dim_and_ADOs_num,
        _check_parity,
        HandleMatrixType,
        _HandleVectorType,
        _HandleSteadyStateMatrix

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
    import LinearAlgebra: I, kron, tr, norm, dot, mul!
    import SparseArrays: sparse, sparsevec, spzeros, SparseVector, SparseMatrixCSC, AbstractSparseMatrix
    import StaticArraysCore: SVector
    import QuantumToolbox:
        _FType,
        QuantumObject,
        Operator,
        SuperOperator,
        TimeDependentOperatorSum,
        _spre,
        _spost,
        spre,
        spost,
        sprepost,
        expect,
        ket2dm,
        liouvillian,
        lindblad_dissipator,
        steadystate,
        ProgressBar,
        next!

    import FastExpm: fastExpm
    import SciMLBase: init, solve, solve!, u_modified!, ODEProblem, FullSpecialize, CallbackSet
    import SciMLOperators: MatrixOperator
    import OrdinaryDiffEqCore: OrdinaryDiffEqAlgorithm
    import OrdinaryDiffEqLowOrderRK: DP5
    import DiffEqCallbacks: PresetTimeCallback, TerminateSteadyState
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
        TimeEvolutionHEOMSol,
        HEOMsolve,
        evolution,   # has been deprecated, throws error only
        SteadyState, # has been deprecated, throws error only
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
    include("steadystate.jl")
    include("power_spectrum.jl")
    include("density_of_states.jl")
end
@reexport using .HeomAPI

end

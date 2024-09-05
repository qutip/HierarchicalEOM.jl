module HierarchicalEOM
import Reexport: @reexport

# sub-module HeomBase for HierarchicalEOM
module HeomBase
    import Pkg
    import LinearAlgebra: BLAS, kron, I
    import SparseArrays: sparse, SparseVector, SparseMatrixCSC
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

# sub-module Bath for HierarchicalEOM
module Bath
    using ..HeomBase
    import Base: show, length, getindex, lastindex, iterate, checkbounds
    import LinearAlgebra: ishermitian, eigvals, I
    import SparseArrays: SparseMatrixCSC
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
    import LinearAlgebra: I, kron, tr
    import SparseArrays: sparse, sparsevec, spzeros, SparseVector, SparseMatrixCSC
    import QuantumToolbox:
        QuantumObject, Operator, SuperOperator, _spre, _spost, spre, spost, sprepost, ket2dm, lindblad_dissipator
    @reexport import QuantumToolbox: expect
    import ProgressMeter: Progress, next!
    import FastExpm: fastExpm

    # for solving time evolution
    import SciMLOperators: MatrixOperator
    import OrdinaryDiffEq: ODEProblem, init, DP5, step!
    import JLD2: jldopen

    # for solving steady state
    import LinearSolve: LinearProblem, init, solve!, UMFPACKFactorization
    import LinearAlgebra: norm
    import OrdinaryDiffEq: solve
    import DiffEqCallbacks: TerminateSteadyState

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
        SteadyState

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
end
@reexport using .HeomAPI

# sub-module Spectrum for HierarchicalEOM
module Spectrum
    using ..HeomBase
    import ..HeomAPI: HEOMSuperOp, ADOs, EVEN, ODD
    import LinearSolve: LinearProblem, init, solve!, UMFPACKFactorization
    import ProgressMeter: Progress, next!

    export PowerSpectrum, DensityOfStates

    include("power_spectrum.jl")
    include("density_of_states.jl")
end
@reexport using .Spectrum

end

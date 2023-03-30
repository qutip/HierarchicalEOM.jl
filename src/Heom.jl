module Heom
    import Reexport: @reexport
    
    export 
        Bath, HeomAPI, Spectrum

    # sub-module HeomBase for Heom
    module HeomBase
        import Pkg
        import LinearAlgebra: BLAS
        import Crayons: Crayon

        include("HeomBase.jl")
    end
    import .HeomBase.versioninfo as versioninfo
    import .HeomBase.print_logo  as print_logo
    
    # sub-module Bath for Heom
    module Bath
        import Base: show, length, getindex, lastindex, iterate, checkbounds
        import LinearAlgebra: I, kron, ishermitian
        import SparseArrays: sparse, SparseMatrixCSC
        import ..HeomBase: isValidMatrixType

        export 
            AbstractBath, BosonBath, FermionBath, Exponent, C,
            AbstractBosonBath, bosonReal, bosonImag, bosonRealImag,
            AbstractFermionBath, fermionAbsorb, fermionEmit,
            spre, spost

        include("Bath.jl")
    end
    @reexport using .Bath

    # sub-module CorrelationFunc for Heom
    module CorrelationFunc
        import ..Bath: BosonBath, FermionBath
        import LinearAlgebra: eigvals

        export 
            Boson_DrudeLorentz_Matsubara, Boson_DrudeLorentz_Pade, 
            Fermion_Lorentz_Matsubara, Fermion_Lorentz_Pade

        include("correlation_functions/CorrelationFunc.jl")
    end
    @reexport using .CorrelationFunc
    
    # sub-module HeomAPI for Heom
    module HeomAPI
        using ..Bath
        import Base: ==, show, length, size, getindex, keys, setindex!, lastindex, iterate, checkbounds, hash, copy
        import Base.Threads: @threads, threadid, nthreads, lock, unlock, SpinLock
        import LinearAlgebra: I, kron, tr
        import SparseArrays: sparse, spzeros, sparsevec, reshape, SparseVector, SparseMatrixCSC, AbstractSparseMatrix
        import ProgressMeter: Progress, next!
        import FastExpm: fastExpm
        import ..HeomBase: PROGBAR_OPTIONS, isValidMatrixType

        # for solving time evolution
        import OrdinaryDiffEq: DiffEqArrayOperator, ODEProblem, init, DP5, step!
        import JLD2: jldopen

        # for solving steady state
        import LinearSolve: LinearProblem, solve, UMFPACKFactorization
        import OrdinaryDiffEq: ODEFunction, SteadyStateProblem, solve, FBDF
        import SteadyStateDiffEq: DynamicSS

        export
            AbstractHEOMMatrix, M_S, M_Boson, M_Fermion, M_Boson_Fermion,
            odd, even,
            ADOs, getRho, getADO, expect,
            Nvec, AbstractHierarchyDict, HierarchyDict, MixHierarchyDict, getIndexEnsemble,
            Propagator, addBosonDissipator, addFermionDissipator, addTerminator,
            evolution, SteadyState

        include("heom_api.jl")
    end
    @reexport using .HeomAPI

    # sub-module Spectrum for Heom
    module Spectrum
        import ..HeomAPI: AbstractHEOMMatrix, ADOs, spre
        import LinearAlgebra: I, kron
        import SparseArrays: sparse, sparsevec
        import LinearSolve: LinearProblem, init, set_A, solve, UMFPACKFactorization
        import ProgressMeter: Progress, next!        
        import ..HeomBase: PROGBAR_OPTIONS, isValidMatrixType

        export spectrum

        include("Spectrum.jl")
    end
    @reexport using .Spectrum

    include("precompile.jl")
end
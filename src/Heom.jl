module Heom
    import Reexport: @reexport
    
    export 
        Bath, HeomBase, Spectrum

    const PROGBAR_OPTIONS = Dict(:barlen=>20, :color=>:green, :showspeed=>true)
    
    # sub-module Bath for Heom
    module Bath
        import LinearAlgebra: eigvals

        export 
            AbstractBath, BosonicBath, FermionicBath,
            pade_NmN

        include("BathBase.jl")
        include("extra_func.jl")
    end
    @reexport using .Bath
    
    # sub-module HeomBase for Heom
    module HeomBase
        import ..Bath: AbstractBath, BosonicBath, FermionicBath
        import Base: size
        import LinearAlgebra: I, kron
        import OrderedCollections: OrderedDict
        import SparseArrays: sparse, reshape, SparseVector, SparseMatrixCSC, AbstractSparseMatrix
        import Distributed: @everywhere, @distributed, procs, nprocs, RemoteChannel, Channel
        import DistributedArrays: distribute, localpart
        import ProgressMeter: Progress, next!
        import ..Heom: PROGBAR_OPTIONS

        export
            AbstractHEOMMatrix, M_Fermion, M_Boson, M_Boson_Fermion, M_CavBath, 
            ADOs, getRho, getADO,
            addDissipator!, spre, spost

        include("HeomBase.jl")
        include("ADOs.jl")
        include("M_fermion.jl")
        include("M_boson.jl")
        include("M_boson_fermion.jl")
        include("M_CavBath.jl")
    end
    @reexport using .HeomBase

    # sub-module evolution for Heom
    module Evolution
        import ..HeomBase: AbstractHEOMMatrix, ADOs
        import OrdinaryDiffEq: ODEProblem, init, DP5, step!
        import SparseArrays: sparse, sparsevec, SparseVector
        import ProgressMeter: Progress, next!
        import ..Heom: PROGBAR_OPTIONS

        export evolution

        include("evolution.jl")
    end
    @reexport using .Evolution

    # sub-module SteadyState for Heom
    module SteadyState
        import ..HeomBase: AbstractHEOMMatrix, ADOs
        import SparseArrays: sparse, sparsevec
        import LinearSolve: LinearProblem, solve, KLUFactorization

        export steadystate

        include("SteadyState.jl")
    end
    @reexport using .SteadyState

    # sub-module Spectrum for Heom
    module Spectrum
        import ..HeomBase: AbstractHEOMMatrix, ADOs, spre
        import LinearAlgebra: I, kron
        import SparseArrays: sparse, sparsevec, SparseVector
        import LinearSolve: LinearProblem, solve, set_A, KLUFactorization
        import ProgressMeter: Progress, next!        
        import ..Heom: PROGBAR_OPTIONS

        export DOS

        include("DensityOfStates.jl")
    end
    @reexport using .Spectrum
end
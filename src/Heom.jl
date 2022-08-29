module Heom
    import Base: size
    import LinearAlgebra: eigvals, I, kron
    import OrderedCollections: OrderedDict
    import OrdinaryDiffEq: ODEProblem, init, DP5, step!
    import SparseArrays: sparse, spzeros, sparsevec, reshape, SparseMatrixCSC, SparseVector, AbstractSparseMatrix
    import ProgressMeter: Progress, next!
    import Distributed: @everywhere, @distributed, procs, nprocs, RemoteChannel, Channel
    import DistributedArrays: distribute, localpart
    
    export AbstractHEOMMatrix, M_fermion, M_boson, M_boson_fermion, M_CavBath, evolution, pade_NmN, Correlation, spre, spost, liouvillian

    PROGBAR_OPTIONS = Dict(:barlen=>20, :color=>:green, :showspeed=>true)

    include("extra_func.jl")
    include("HeomBase.jl")
    include("M_fermion.jl")
    include("M_boson.jl")
    include("M_boson_fermion.jl")
    include("M_CavBath.jl")
end
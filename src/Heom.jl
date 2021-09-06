module Heom
    import Base: size
    import LinearAlgebra: eigvals, I, kron
    import DifferentialEquations: ODEProblem, init, Tsit5, step!
    import SparseArrays: sparse, spzeros, SparseMatrixCSC, SparseVector, AbstractSparseMatrix
    import QuantumOptics: AbstractOperator, dagger
    import ProgressMeter: Progress, next!

    export AbstractHEOMMatrix, M_fermion, M_boson, M_boson_fermion, evolution, pade_NmN, Correlation, spre, spost, liouvillian

    PROGBAR_OPTIONS = Dict(:barlen=>20, :color=>:green, :showspeed=>true)

    include("extra_func.jl")
    include("HeomBase.jl")
    include("M_fermion.jl")
    include("M_boson.jl")
    include("M_boson_fermion.jl")
end
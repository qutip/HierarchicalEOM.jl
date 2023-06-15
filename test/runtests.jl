using HierarchicalEOM
using Test
using SparseArrays
using LinearAlgebra

include("utils.jl")

@testset "Print version information" begin
    @test HierarchicalEOM.versioninfo() == nothing
end

@testset "Bath and Exponent" begin
    include("bath.jl")
end

@testset "Correlation functions" begin
    include("corr_func.jl")
end

include("heom_api.jl")

include("phys_analysis.jl")
using QuantumOptics
using HierarchicalEOM
using Test
using SparseArrays
using LinearAlgebra

const GROUP = get(ENV, "GROUP", "All")
const HAS_EXTENSIONS = isdefined(Base, :get_extension)

include("utils.jl")

if GROUP == "All" || GROUP == "Core"

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
end

if (GROUP == "All" || GROUP == "HierarchicalEOM_QOExt") && HAS_EXTENSIONS
    @testset "QuantumOptics Extension" begin
        include("QOExt.jl")
    end
end
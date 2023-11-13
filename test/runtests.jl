using Test, Pkg
using HierarchicalEOM
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

    @testset "Bath correlation functions" begin
        include("bath_corr_func.jl")
    end

    include("heom_api.jl")

    include("phys_analysis.jl")
end

if GROUP == "HierarchicalEOM_CUDAExt"
    Pkg.activate("CUDA")
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    Pkg.instantiate()
    @testset "CUDA Extension" begin
        include("CUDA/CUDAExt.jl")
    end
end

if (GROUP == "All" || GROUP == "HierarchicalEOM_QOExt") && HAS_EXTENSIONS
    Pkg.activate("QuantumOptics")
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    Pkg.instantiate()
    @testset "QuantumOptics Extension" begin
        include("QuantumOptics/QOExt.jl")
    end
end
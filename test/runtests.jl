using HierarchicalEOM
using Test, Pkg, SparseArrays, LinearAlgebra

const GROUP = get(ENV, "GROUP", "All")

include("utils.jl")

if (GROUP == "All") || GROUP == "Core"

    @testset "Print version information" begin
        @test HierarchicalEOM.versioninfo() == nothing
    end

    @testset "Bath and Exponent" begin
        include("bath.jl")
    end

    @testset "Bath correlation functions" begin
        include("bath_corr_func.jl")
    end

    @testset "Auxiliary density operators" begin
        include("ADOs.jl")
    end

    @testset "HEOM superoperator" begin
        include("HEOMSuperOp.jl")
    end

    @testset "M_S" begin
        include("M_S.jl")
    end

    @testset "M_Boson" begin
        include("M_Boson.jl")
    end

    @testset "M_Boson (RWA)" begin
        include("M_Boson_RWA.jl")
    end

    @testset "M_Fermion" begin
        include("M_Fermion.jl")
    end

    @testset "M_Boson_Fermion" begin
        include("M_Boson_Fermion.jl")
    end

    @testset "Hierarchy Dictionary" begin
        include("hierarchy_dictionary.jl")
    end

    @testset "Stationary state" begin
        include("stationary_state.jl")
    end

    @testset "Time evolution" begin
        include("time_evolution.jl")
    end

    @testset "Power spectrum" begin
        include("power_spectrum.jl")
    end

    @testset "Density of states" begin
        include("density_of_states.jl")
    end
end

if GROUP == "HierarchicalEOM_CUDAExt"
    Pkg.activate("CUDA")
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    Pkg.instantiate()
    @testset "CUDA Extension" begin
        include("CUDA/CUDAExt.jl")
    end
end

if (GROUP == "All") || (GROUP == "HierarchicalEOM_QOExt")
    Pkg.activate("QuantumOptics")
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    Pkg.instantiate()
    @testset "QuantumOptics Extension" begin
        include("QuantumOptics/QOExt.jl")
    end
end

if (GROUP == "All") || (GROUP == "HierarchicalEOM_QuantumToolboxExt")
    Pkg.activate("QuantumToolbox")
    Pkg.add(url="https://github.com/albertomercurio/QuantumToolbox.jl")
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    Pkg.instantiate()
    @testset "QuantumToolbox Extension" begin
        include("QuantumToolbox/QTExt.jl")
    end
end
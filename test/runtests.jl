using Pkg
using Test
using SafeTestsets

const GROUP = get(ENV, "GROUP", "All")

if (GROUP == "All") || GROUP == "Core"

    @testset "Print version information" begin
        using HierarchicalEOM
        @test HierarchicalEOM.versioninfo() == nothing
    end

    @time @safetestset "Bath and Exponent" include("bath.jl")
    @time @safetestset "Bath correlation functions" include("bath_corr_func.jl")
    @time @safetestset "Auxiliary density operators" include("ADOs.jl")
    @time @safetestset "HEOM superoperator" include("HEOMSuperOp.jl")
    @time @safetestset "M_S" include("M_S.jl")
    @time @safetestset "M_Boson" include("M_Boson.jl")
    @time @safetestset "M_Boson (RWA)" include("M_Boson_RWA.jl")
    @time @safetestset "M_Fermion" include("M_Fermion.jl")
    @time @safetestset "M_Boson_Fermion" include("M_Boson_Fermion.jl")
    @time @safetestset "Hierarchy Dictionary" include("hierarchy_dictionary.jl")
    @time @safetestset "Stationary state" include("stationary_state.jl")
    @time @safetestset "Time evolution" include("time_evolution.jl")
    @time @safetestset "Power spectrum" include("power_spectrum.jl")
    @time @safetestset "Density of states" include("density_of_states.jl")
end

if GROUP == "HierarchicalEOM_CUDAExt"
    Pkg.activate("CUDA")
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    Pkg.instantiate()
    @time @safetestset "CUDA Extension" include("CUDA/CUDAExt.jl")
end

if (GROUP == "All") || (GROUP == "HierarchicalEOM_QuantumOpticsExt")
    Pkg.activate("QuantumOptics")
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    Pkg.instantiate()
    @time @safetestset "QuantumOptics Extension" include("QuantumOptics/QOExt.jl")
end

if (GROUP == "All") || (GROUP == "HierarchicalEOM_QuantumToolboxExt")
    Pkg.activate("QuantumToolbox")
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    Pkg.instantiate()
    @time @safetestset "QuantumToolbox Extension" include("QuantumToolbox/QTExt.jl")
end
using Pkg
using Test
using HierarchicalEOM
using LinearAlgebra, SparseArrays

const GROUP = get(ENV, "GROUP", "All")

include("test_utils.jl")

# Put Core tests in alphabetical order
core_tests = [
    "ADOs.jl",
    "bath.jl",
    "bath_corr_func.jl",
    "density_of_states.jl",
    "HEOMSuperOp.jl",
    "hierarchy_dictionary.jl",
    "M_Boson.jl",
    "M_Boson_Fermion.jl",
    "M_Boson_RWA.jl",
    "M_Fermion.jl",
    "M_S.jl",
    "power_spectrum.jl",
    "stationary_state.jl",
    "time_evolution.jl"
]

HierarchicalEOM.versioninfo()

if (GROUP == "All") || (GROUP == "Core")
    for test in core_tests
        include(test)
    end
end

if GROUP == "HierarchicalEOM_CUDAExt"
    Pkg.activate("CUDA")
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    Pkg.instantiate()
    include("CUDA/CUDAExt.jl")
end

if (GROUP == "All") || (GROUP == "HierarchicalEOM_QuantumOpticsExt")
    Pkg.activate("QuantumOptics")
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    Pkg.instantiate()
    include("QuantumOptics/QOExt.jl")
end

if (GROUP == "All") || (GROUP == "HierarchicalEOM_QuantumToolboxExt")
    Pkg.activate("QuantumToolbox")
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    Pkg.instantiate()
    include("QuantumToolbox/QTExt.jl")
end
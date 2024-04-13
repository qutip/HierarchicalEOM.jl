using Pkg
using Test
using HierarchicalEOM
using LinearAlgebra, SparseArrays

const GROUP = get(ENV, "GROUP", "All")

const testdir = dirname(@__FILE__)

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
        include(joinpath(testdir, test))
    end
end

if (GROUP == "All") || (GROUP == "Code_Quality")
    Pkg.add(["Aqua", "JET"])
    include(joinpath(testdir, "aqua.jl"))
    include(joinpath(testdir, "jet.jl"))
end

if GROUP == "HierarchicalEOM_CUDAExt"
    Pkg.add("CUDA")
    include(joinpath(testdir, "CUDA_Ext.jl"))
end

if (GROUP == "All") || (GROUP == "QuantumOptics_Ext")
    Pkg.add("QuantumOptics")
    include(joinpath(testdir, "QuantumOpticsExt.jl"))
end

if (GROUP == "All") || (GROUP == "QuantumToolbox_Ext")
    Pkg.add("QuantumToolbox")
    include(joinpath(testdir, "QuantumToolboxExt.jl"))
end
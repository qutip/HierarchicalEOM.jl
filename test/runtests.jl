using Pkg
using Test
using HierarchicalEOM
using HierarchicalEOM: spre, spost
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
    "time_evolution.jl",
]

if (GROUP == "All") || (GROUP == "Code_Quality")
    Pkg.add(["Aqua", "JET"])

    HierarchicalEOM.versioninfo()
    include(joinpath(testdir, "aqua.jl"))
    include(joinpath(testdir, "jet.jl"))
end

if (GROUP == "All") || (GROUP == "Core")
    GROUP == "All" ? nothing : HierarchicalEOM.versioninfo()
    for test in core_tests
        include(joinpath(testdir, test))
    end
end

if GROUP == "CUDA_Ext"
    Pkg.add("CUDA")

    HierarchicalEOM.versioninfo()
    include(joinpath(testdir, "CUDAExt.jl"))
end

if (GROUP == "All") || (GROUP == "QuantumOptics_Ext")
    Pkg.add("QuantumOptics")

    GROUP == "All" ? nothing : HierarchicalEOM.versioninfo()
    include(joinpath(testdir, "QuantumOpticsExt.jl"))
end

if (GROUP == "All") || (GROUP == "QuantumToolbox_Ext")
    Pkg.add("QuantumToolbox")

    GROUP == "All" ? nothing : HierarchicalEOM.versioninfo()
    include(joinpath(testdir, "QuantumToolboxExt.jl"))
end

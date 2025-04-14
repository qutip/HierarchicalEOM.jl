using Pkg
using Test

# Importing only the necessary functions to keep track the re-export of the functions
import SparseArrays: spzeros, nnz

const GROUP = get(ENV, "GROUP", "All")

const testdir = dirname(@__FILE__)

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

if (GROUP == "All") || (GROUP == "Core")
    using HierarchicalEOM
    import JLD2: jldopen

    HierarchicalEOM.about()

    include(joinpath(testdir, "test_utils.jl"))

    for test in core_tests
        include(joinpath(testdir, test))
    end
end

if (GROUP == "All") || (GROUP == "Code-Quality")
    Pkg.activate("code-quality")
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    Pkg.instantiate()

    using HierarchicalEOM
    using Aqua, JET

    include(joinpath(testdir, "code-quality", "code_quality.jl"))
end

if (GROUP == "CUDA_Ext")# || (GROUP == "All")
    Pkg.activate("gpu")
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    Pkg.instantiate()

    using HierarchicalEOM
    using CUDA, LinearSolve
    CUDA.allowscalar(false) # Avoid unexpected scalar indexing

    HierarchicalEOM.about()
    CUDA.versioninfo()

    include(joinpath(testdir, "test_utils.jl"))
    include(joinpath(testdir, "gpu", "CUDAExt.jl"))
end

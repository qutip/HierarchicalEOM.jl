using Test
using TestItemRunner
using Pkg

const GROUP_LIST = String["All", "Core", "Code-Quality", "CUDA_Ext"]

const GROUP = get(ENV, "GROUP", "All")
(GROUP in GROUP_LIST) || throw(ArgumentError("Unknown GROUP = $GROUP"))

if (GROUP == "All") || (GROUP == "Core")
    import HierarchicalEOM

    HierarchicalEOM.about()

    println("\nStart running Core tests...\n")
    @run_package_tests verbose=true
end

########################################################################
# Use traditional Test.jl instead of TestItemRunner.jl for other tests #
########################################################################

const testdir = dirname(@__FILE__)

if (GROUP == "All") || (GROUP == "Code-Quality")
    Pkg.activate("code-quality")
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    Pkg.instantiate()

    using HierarchicalEOM
    using Aqua, JET

    (GROUP == "Code-Quality") && HierarchicalEOM.about() # print version info. for code quality CI in GitHub

    include(joinpath(testdir, "code-quality", "code_quality.jl"))
end

if (GROUP == "CUDA_Ext")# || (GROUP == "All")
    Pkg.activate("gpu")
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    Pkg.instantiate()

    using HierarchicalEOM
    using LinearSolve
    using CUDA
    CUDA.allowscalar(false) # Avoid unexpected scalar indexing

    HierarchicalEOM.about()
    CUDA.versioninfo()

    include(joinpath(testdir, "gpu", "CUDAExt.jl"))
end

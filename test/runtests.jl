using ParallelTestRunner
using Pkg

include("group_list.jl")

const GROUP = get(ENV, "GROUP", "All")
(GROUP in GROUP_LIST) || throw(ArgumentError("Unknown GROUP = $GROUP\nAvailable test GROUP are:\n$SHOW_GROUP_LIST\n"))

# function to set up the environment for subtests
function setup_subtest_env(path::String)
    Pkg.activate(path)

    # `cd` into `path` first; otherwise `Pkg.develop` writes an absolute (machine-specific) path to Project.toml
    cd(path) do
        Pkg.develop(PackageSpec(path = relpath(dirname(@__DIR__), path))) # must `develop` otherwise the code coverage for libraries (like CUDA) will not be collected
    end
    Pkg.update()
    return nothing
end

######################
# Main package tests #
######################
if (GROUP == "All") || (GROUP == "Main")
    testsuite = find_tests(joinpath(testdir, "main-test"))

    import HierarchicalEOM
    HierarchicalEOM.about()
    println("[Tests for GROUP = $GROUP]")
    runtests(HierarchicalEOM, ARGS; testsuite)
end

######################
# Code Quality tests #
######################

if (GROUP == "All") || (GROUP == "Code-Quality")
    path = joinpath(testdir, "code-quality")
    setup_subtest_env(path)

    using HierarchicalEOM
    using Aqua, JET

    (GROUP == "Code-Quality") && HierarchicalEOM.about() # print version info. for code quality CI in GitHub

    println("[Tests for GROUP = $GROUP]")
    include(joinpath(path, "code_quality.jl"))
end

###################
# Extension tests #
###################
if GROUP ∈ EXTENSION_LIST
    path = EXTENSION_PATH[GROUP]
    setup_subtest_env(path)

    println("[Tests for GROUP = $GROUP]")
    include(joinpath(path, "runtests.jl"))
end

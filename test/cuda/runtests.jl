using HierarchicalEOM
using CUDA
using ParallelTestRunner

HierarchicalEOM.about()
CUDA.versioninfo()

testsuite = find_tests(dirname(@__FILE__))

runtests(HierarchicalEOM, ARGS; testsuite)

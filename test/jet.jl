import JET

@testset "Code quality (JET.jl)" begin

JET.test_package(HierarchicalEOM; target_defined_modules=true, ignore_missing_comparison=true)

end
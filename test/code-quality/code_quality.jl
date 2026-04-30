@testset "Code quality" verbose = true begin
    @testset "Aqua.jl" begin
        Aqua.test_all(HierarchicalEOM; ambiguities = false, unbound_args = false)
    end

    @testset "JET.jl" begin
        JET.test_package(HierarchicalEOM; target_modules = (HierarchicalEOM,), ignore_missing_comparison = true)
    end
end

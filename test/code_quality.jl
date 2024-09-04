import Aqua
import JET

@testset "Code quality" verbose = true begin
    @testset "Aqua.jl" begin
        Aqua.test_all(HierarchicalEOM; ambiguities = false)
    end

    @testset "JET.jl" begin
        JET.test_package(HierarchicalEOM; target_defined_modules = true, ignore_missing_comparison = true)
    end
end
import Aqua

@testset "Code quality (Aqua.jl)" begin
    Aqua.test_all(HierarchicalEOM; ambiguities = false)
end

# System Hamiltonian and initial state
Hsys = 0.25 * [1 0; 0 -1] + 0.5 * [0 1; 1 0]
ρ0   = [1 0; 0 0]

# Bath properties:
λ = 0.1
W = 0.5
T = 0.5
N = 2
Q = [1 0; 0 -1]  # System-bath coupling operator
bath = Boson_DrudeLorentz_Pade(Q, λ, W, T, N)

# Heom liouvillian superoperator matrix
tier = 5
L = M_Boson(Hsys, tier, bath; verbose=false)

ρs = getRho(SteadyState(L, ρ0; verbose=false))
@testset "Steady state" begin
    ρ1 = getRho(SteadyState(L; verbose=false))
    @test _is_Matrix_approx(ρ1, ρs)
end

@testset "Time evolution" begin
    Δt    = 10
    steps = 10
    tlist = 0:Δt:(Δt * steps)
    ρ_list_p = getRho.(evolution(L, ρ0, Δt, steps; verbose=false))  # using the method based on propagator
    ρ_list_e = getRho.(evolution(L, ρ0, tlist; verbose=false))      # using the method based on ODE solver
    for i in 1:(steps + 1)
        @test _is_Matrix_approx(ρ_list_p[i], ρ_list_e[i])
    end
    @test _is_Matrix_approx(ρs, ρ_list_p[end])
    @test _is_Matrix_approx(ρs, ρ_list_e[end])
end
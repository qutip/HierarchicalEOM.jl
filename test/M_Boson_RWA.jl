@time @testset "M_Boson (RWA)" begin

# Test Boson-type HEOM Liouvillian superoperator matrix under rotating wave approximation
ωq = 1.1
Λ  = 0.01
Γ  = 0.02
tier = 3
Hsys_rwa = 0.5 * ωq * [1 0; 0 -1]
op_rwa   = [0 0; 1 0]
ρ0       = 0.5 * [1 1; 1 1]

tlist = 0:1:20
d     = 1im * √(Λ * (2 * Γ - Λ)) # non-Markov regime

B_rwa = BosonBathRWA(op_rwa, [0], [Λ - 1im * ωq], [0.5 * Γ * Λ], [Λ + 1im * ωq])
L = M_Boson(Hsys_rwa, tier, B_rwa; verbose=false)
ados_list = evolution(L, ρ0, tlist; reltol=1e-10, abstol=1e-12, verbose=false)

for (i, t) in enumerate(tlist)
    ρ_rwa = getRho(ados_list[i])

    # analytical result
    Gt = exp(-1 * (Λ / 2 + 1im * ωq) * t) * (cosh(t * d / 2) + Λ * sinh(t * d / 2) / d)
    
    @test ρ_rwa[1, 1] ≈ abs(Gt) ^ 2 * ρ0[1, 1]
    @test ρ_rwa[1, 2] ≈ Gt * ρ0[1, 2]
end
end
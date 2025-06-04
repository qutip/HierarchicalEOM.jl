@testitem "Hierarchy Dictionary" begin
    using SparseArrays

    # calculate current for a given ADOs
    # bathIdx: 1 means 1st fermion bath (bath_L); 2 means 2nd fermion bath (bath_R)
    function Ic(ados, M::M_Boson_Fermion, bathIdx::Int)
        # the hierarchy dictionary
        HDict = M.hierarchy

        # we need all the indices of ADOs for the first level
        idx_list = HDict.Flvl2idx[1]
        I = 0.0im
        for idx in idx_list
            ρ1 = ados[idx]  # 1st-level ADO

            # find the corresponding bath index (α) and exponent term index (k)
            nvec_b, nvec_f = HDict.idx2nvec[idx]
            if nvec_b.level == 0
                for (α, k, _) in getIndexEnsemble(nvec_f, HDict.fermionPtr)
                    if α == bathIdx
                        exponent = M.Fbath[α][k]
                        if exponent.types == "fA"     # fermion-absorption
                            I += tr(exponent.op' * ρ1)
                        elseif exponent.types == "fE" # fermion-emission
                            I -= tr(exponent.op' * ρ1)
                        end
                        break
                    end
                end
            end
        end

        eV_to_Joule = 1.60218E-19  # unit conversion

        # (e / ħ) * I  [change unit to μA] 
        return 1.519270639695384E15 * real(1im * I) * eV_to_Joule * 1E6
    end

    Btier = 2
    Ftier = 2
    Nb = 3
    Nf = 3
    threshold = 1e-5

    Γ = 0.0025
    Dα = 30
    Λ = 0.0025
    ωcα = 0.2
    μL = 0.5
    μR = -0.5
    kT = 0.025

    Hsys = Qobj([
        0 0 0 0
        0 0.2 0 0
        0 0 0.208 0.04
        0 0 0.04 0.408
    ])

    cop = Qobj([
        0 1 0 0
        1 0 0 0
        0 0 0 1
        0 0 1 0
    ])

    dop = Qobj([
        0 0 1 0
        0 0 0 1
        0 0 0 0
        0 0 0 0
    ])

    bbath = Boson_DrudeLorentz_Matsubara(cop, Λ, ωcα, kT, Nb)
    fbath = [Fermion_Lorentz_Pade(dop, Γ, μL, Dα, kT, Nf), Fermion_Lorentz_Pade(dop, Γ, μR, Dα, kT, Nf)]

    L = M_Boson_Fermion(Hsys, Btier, Ftier, bbath, fbath; threshold = threshold, verbose = false)
    L = addTerminator(L, bbath)

    @test size(L) == (1696, 1696)
    @test L.N == 106
    @test nnz(L.data.A) == nnz(L(0)) == 13536

    ados = steadystate(L; verbose = false)
    @test Ic(ados, L, 1) ≈ 0.2883390125832726
    nvec_b, nvec_f = L.hierarchy.idx2nvec[1]
    @test_throws ErrorException getIndexEnsemble(nvec_f, L.hierarchy.bosonPtr)
    @test_throws ErrorException getIndexEnsemble(nvec_b, L.hierarchy.fermionPtr)
end

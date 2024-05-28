function _is_Matrix_approx(M1, M2; atol = 1.0e-6)
    s1 = size(M1)
    s2 = size(M2)
    if s1 == s2
        m, n = s1
        for i in 1:m
            for j in 1:n
                if !isapprox(abs(M1[i, j]), abs(M2[i, j]), atol = atol)
                    return false
                end
            end
        end
        return true
    else
        return false
    end
end

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

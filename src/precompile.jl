import PrecompileTools

PrecompileTools.@setup_workload begin
    # Putting some things in `setup` can reduce the size of the
    # precompile file and potentially make loading faster.
    d  = [0 1; 0 0]
    Hs = d' * d
    ρ0 = [1 0; 0 0]
    bB = Boson_DrudeLorentz_Pade(d' * d, 1, 1., 0.05, 3)
    fB = Fermion_Lorentz_Pade(d, 1., 0, 1., 0.05, 3)

    PrecompileTools.@compile_workload begin
        # all calls in this block will be precompiled, regardless of whether
        # they belong to your package or not (on Julia 1.8 and higher)
        
        # precompile bath and exponents
        bB[1:end]
        fB[1:end]
        for b in bB nothing end
        for b in fB nothing end

        # precompile HEOM matrices
        Ms   = M_S(Hs; verbose=false)
        Mb   = M_Boson(Hs, 2, bB; verbose=false, threshold=1e-6)
        Mfe  = M_Fermion(Hs, 5, fB, EVEN; verbose=false, threshold=1e-8)
        Mfo  = M_Fermion(Hs, 5, fB, ODD;  verbose=false, threshold=1e-8)
        Mbfe = M_Boson_Fermion(Hs, 2, 2, bB, fB; verbose=false, threshold=1e-6)

        # precompile Steadystate
        ados1 = SteadyState(Mfe; verbose=false)
        ados1 = SteadyState(Mfe, ρ0; verbose=false)
        E1 = Expect(Hs, ados1)

        # precompile evolution
        ados2  = evolution(Mb, ρ0, 0:1:1; verbose=false)
        ados2  = evolution(Mb, ρ0, 1, 1; verbose=false)
        E2 = Expect(Hs, ados2)

        # precompile ADOs for the support of Base functions 
        ados3 = ADOs(zeros(8), 2)
        length(ados3)
        getRho(ados3)
        getADO(ados3, 1)
        ados3[end]
        ados3[:]
        for ad in ados3 nothing end

        # precompile Spectrum functions
        psd = PowerSpectrum(  Mfo, ados1, d, [1]; verbose=false)
        dos = DensityOfStates(Mfo, ados1, d, [1]; verbose=false)
    end
end
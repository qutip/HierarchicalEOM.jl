import PrecompileTools

PrecompileTools.@setup_workload begin
    # Putting some things in `setup` can reduce the size of the
    # precompile file and potentially make loading faster.
    op = [0.1 0.2; 0.2 0.4]
    bB = Boson_DrudeLorentz_Pade(op, 1, 1., 1., 3)
    fB = Fermion_Lorentz_Pade(op, 1., 1., 1., 1., 2)

    PrecompileTools.@compile_workload begin
        # all calls in this block will be precompiled, regardless of whether
        # they belong to your package or not (on Julia 1.8 and higher)
        
        # precompile bath and exponents
        @info "Precompiling Bath and Exponent..."
        bB[1:end]
        fB[1:end]
        for b in bB nothing end
        for b in fB nothing end

        # precompile HEOM matrices
        @info "Precompiling HEOM Liouvillian superoperator matrices..."
        Ms   = M_S(op; verbose=false)
        Mb   = M_Boson(op, 2, bB; verbose=false, threshold=1e-1)
        Mfo  = M_Fermion(op, 2, fB, :odd; verbose=false, threshold=1e-1)
        Mbfe = M_Boson_Fermion(op, 2, 2, bB, fB; verbose=false, threshold=1e-1)

        # precompile Steadystate
        @info "Precompiling steady state solver..."
        ados1 = SteadyState(Mb; verbose=false)
        ados1 = SteadyState(Mb, [1. 0.; 0. 0.]; verbose=false)
        E1 = Expect(op, ados1)

        # precompile evolution
        @info "Precompiling time evolution solver..."
        ados2  = evolution(Mb, [1. 0.; 0. 0.], 0:1:1; verbose=false)
        ados2  = evolution(Mb, [1. 0.; 0. 0.], 1, 1; verbose=false)
        E2 = Expect(op, ados2)

        # precompile ADOs for the support of Base functions 
        @info "Precompiling Auxiliary Density Operators (ADOs)..."
        ados3 = ADOs(zeros(8),  2)
        length(ados3)
        getRho(ados3)
        getADO(ados3, 1)
        ados3[end]
        ados3[:]
        for ad in ados3 nothing end

        # precompile Spectrum functions
        @info "Precompiling solvers for calculating spectrum..."
        psd = spectrum(Mb,  [1. 0.; 0. 0.], op, [1]; verbose=false)
        dos = spectrum(Mfo, [1. 0.; 0. 0.], op, [1]; verbose=false)
    end
    @info "HEOM precompilation complete"
end
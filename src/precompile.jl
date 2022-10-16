import SnoopPrecompile: @precompile_setup, @precompile_all_calls

@precompile_setup begin
    # Putting some things in `setup` can reduce the size of the
    # precompile file and potentially make loading faster.
    op = [0.1 0.2; 0.3 0.4]
    bB = Boson_DrudeLorentz_Pade(op, 1, 1., 1., 3)
    fB = Fermion_Lorentz_Pade(op, 1., 1., 1., 1., 2)

    @precompile_all_calls begin
        # all calls in this block will be precompiled, regardless of whether
        # they belong to your package or not (on Julia 1.8 and higher)
        
        # precompile bath and exponents
        @info "Precompiling Bath and Exponent..."
        bB[1:end]
        fB[1:end]
        for b in bB nothing end
        for b in fB nothing end

        # precompile Heom matrices
        @info "Precompiling Heom liouvillian superoperator matrices..."
        Mb   = M_Boson(op, 2, bB, verbose=false)
        Mfo  = M_Fermion(op, 2, fB, :odd, verbose=false)
        Mbfe = M_Boson_Fermion(op, 2, 2, bB, fB, verbose=false)

        # precompile Steadystate
        @info "Precompiling steady state solver..."
        SteadyState(Mb; verbose=false)
        SteadyState(Mb, [1. 0.; 0. 0.]; verbose=false)

        # precompile evolution
        @info "Precompiling time evolution solver..."
        adosb  = evolution(Mb, [1. 0.; 0. 0.], 0:1:1; verbose=false)
        adosb  = evolution(Mb, [1. 0.; 0. 0.], 1, 1; verbose=false)

        # precompile ADOs for the support of Base functions 
        @info "Precompiling Auxiliary Density Operators (ADOs)..."
        ados1 = ADOs(zeros(8); Nb = 2)
        ados2 = ADOs(zeros(16), Nb = 2, Nf = 2)
        length(ados1)
        length(ados2)
        getRho(ados1)
        getADO(ados1, 1)
        getADO(ados2, 1, 1)
        ados1[end]
        ados1[:]
        ados2[1, :]
        ados2[:, 1]
        ados2[end, end]
        for ad in ados1 nothing end

        # precompile Spectrum functions
        @info "Precompiling solvers for calculating spectrum..."
        PSD(Mb,  [1. 0.; 0. 0.], op, [1]; verbose=false)
        DOS(Mfo, [1. 0.; 0. 0.], op, [1]; verbose=false)
    end
    GC.gc() # clean garbage collection
    @info "Heom precompilation complete"
end
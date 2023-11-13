@doc raw"""
    DensityOfStates(M, ρ, d_op, ωlist; solver, verbose, filename, SOLVEROptions...)
Calculate density of states for the fermionic system in frequency domain.

```math
    \pi A(\omega)=\textrm{Re}\left\{\int_0^\infty dt \left[\langle d(t) d^\dagger(0)\rangle^* + \langle d^\dagger(t) d(0)\rangle \right] e^{-i\omega t}\right\},
```

# Parameters
- `M::AbstractHEOMLSMatrix` : the HEOMLS matrix which acts on `ODD`-parity operators.
- `ρ` :  the system density matrix or the auxiliary density operators.
- `d_op` : The annihilation operator (``d`` as shown above) acting on the fermionic system.
- `ωlist::AbstractVector` : the specific frequency points to solve.
- `solver` : solver in package `LinearSolve.jl`. Default to `UMFPACKFactorization()`.
- `verbose::Bool` : To display verbose output and progress bar during the process or not. Defaults to `true`.
- `filename::String` : If filename was specified, the value of spectrum for each ω will be saved into the file "filename.txt" during the solving process.
- `SOLVEROptions` : extra options for solver 

For more details about solvers and extra options, please refer to [`LinearSolve.jl`](http://linearsolve.sciml.ai/stable/)

# Returns
- `dos::AbstractVector` : the list of density of states corresponds to the specified `ωlist`
"""
@noinline function DensityOfStates(
        M::AbstractHEOMLSMatrix, 
        ρ,
        d_op,
        ωlist::AbstractVector; 
        solver=UMFPACKFactorization(), 
        verbose::Bool = true,
        filename::String = "",
        SOLVEROptions...
    )

    Size = size(M, 1)

    # check M
    if M.parity == EVEN
        error("The HEOMLS matrix M must be acting on `ODD`-parity operators.")
    end

    # Handle ρ
    if typeof(ρ) == ADOs  # ρ::ADOs
        ados = ρ
        if ados.parity != EVEN
            error("The parity of ρ must be `EVEN`.")
        end
    else
        ados = ADOs(ρ, M.N)
    end
    _check_sys_dim_and_ADOs_num(M, ados)

    # Handle d_op
    _tr = Tr(M.dim, M.N)
    d_normal = HEOMSuperOp(d_op,  M, ODD)
    d_dagger = HEOMSuperOp(d_op', M, ODD)
    b_m = _HandleVectorType(typeof(M.data), (d_normal * ados).data)
    b_p = _HandleVectorType(typeof(M.data), (d_dagger * ados).data)
    _tr_d_normal = _tr * d_normal.data
    _tr_d_dagger = _tr * d_dagger.data

    SAVE::Bool = (filename != "")
    if SAVE
        FILENAME = filename * ".txt"
        if isfile(FILENAME)
            error("FILE: $(FILENAME) already exist.")
        end
    end
    
    ElType = eltype(M)
    ωList  = _HandleFloatType(ElType, ωlist)
    Length = length(ωList)
    Aω = Vector{Float64}(undef, Length)

    if verbose
        print("Calculating density of states in frequency domain...\n")
        flush(stdout)
        prog = Progress(Length; desc="Progress : ", PROGBAR_OPTIONS...)
    end
    i  = convert(ElType, 1im)
    I_total = _HandleIdentityType(typeof(M.data), Size)
    Iω      = i * ωList[1] * I_total
    cache_m = init(LinearProblem(M.data - Iω, b_m),  solver, SOLVEROptions...)
    cache_p = init(LinearProblem(M.data + Iω, b_p), solver, SOLVEROptions...)
    sol_m   = solve!(cache_m)
    sol_p   = solve!(cache_p)
    @inbounds for (j, ω) in enumerate(ωList)
        if j > 1
            Iω = i * ω * I_total

            cache_m.A = M.data - Iω
            sol_m = solve!(cache_m)
            
            cache_p.A = M.data + Iω
            sol_p = solve!(cache_p)
        end
        
        # trace over the Hilbert space of system (expectation value)
        Aω[j] = -1 * (
            real(_tr_d_normal * _HandleVectorType(sol_p.u, false)) + 
            real(_tr_d_dagger * _HandleVectorType(sol_m.u, false))
        )

        if SAVE
            open(FILENAME, "a") do file
                write(file, "$(Aω[j]),\n")
            end
        end

        if verbose
            next!(prog)
        end 
    end
    if verbose
        println("[DONE]")
    end

    return Aω
end
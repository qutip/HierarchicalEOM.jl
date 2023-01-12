"""
    spectrum(M, ρ, op, ω_list; solver, verbose, filename, SOLVEROptions...)
Calculate spectrum for the system.

# To calculate spectrum for bosonic systems (usually known as power spectrum):
```math
\\pi S(\\omega)=\\textrm{Re}\\left\\{\\int_0^\\infty dt \\langle a^\\dagger(t) a(0)\\rangle e^{-i\\omega t}\\right\\},
```
remember to set the parameters: 
- `M::AbstractHEOMMatrix`: should be either `:none` or `:even` parity
- `op`: the (annihilation) operator ``a`` for bosonic system as shown above 

# To calculate spectrum for fermionic systems (usually known as density of states):
```math
    \\pi A(\\omega)=\\textrm{Re}\\left\\{\\int_0^\\infty dt \\left[\\langle d(t) d^\\dagger(0)\\rangle^* + \\langle d^\\dagger(t) d(0)\\rangle \\right] e^{-i\\omega t}\\right\\},
```
remember to set the parameters: 
- `M::AbstractHEOMMatrix`: should be `:odd` parity
- `op`: the (annihilation) operator ``d`` for fermionic system as shown above 

# Parameters
- `M::AbstractHEOMMatrix` : the matrix given from HEOM model.
- `ρ` :  the system density matrix or the auxiliary density operators.
- `op` : The annihilation operator acting on the system.
- `ω_list::AbstractVector` : the specific frequency points to solve.
- `solver` : solver in package `LinearSolve.jl`. Default to `UMFPACKFactorization()`.
- `verbose::Bool` : To display verbose output and progress bar during the process or not. Defaults to `true`.
- `filename::String` : If filename was specified, the value of spectrum for each ω will be saved into the file during the solving process.
- `SOLVEROptions` : extra options for solver 

For more details about solvers and extra options, please refer to [`LinearSolve.jl`](http://linearsolve.sciml.ai/stable/)

# Returns
- `spec::AbstractVector` : the spectrum list corresponds to the specified `ω_list`
"""
function spectrum(
        M::AbstractHEOMMatrix, 
        ρ, 
        op, 
        ω_list::AbstractVector; 
        solver=UMFPACKFactorization(), 
        verbose::Bool = true,
        filename::String = "",
        SOLVEROptions...
    )
    if (filename != "") && isfile(filename)
        error("FILE: $(filename) already exist.")
    end

    Size, = size(M)

    # check ρ
    if typeof(ρ) == ADOs  # ρ::ADOs
        if (M.dim != ρ.dim)
            error("The system dimension between M and ρ are not consistent.")
        end
        if (M.N != ρ.N)
            error("The number N between M and ρ are not consistent.")
        end

        b = ρ.data
        
    elseif isValidMatrixType(ρ, M.dim)
        v = sparsevec(ρ)
        b = sparsevec(v.nzind, v.nzval, Size)
    else
        error("Invalid matrix \"ρ\".")
    end

    # check dimension of op
    if !isValidMatrixType(op, M.dim)
        error("Invalid matrix \"op\".")
    end

    # check parity and calculate spectrum
    if (M.parity == :odd)
        return _density_of_states(M, b, op, ω_list; 
            solver = solver, 
            verbose = verbose,
            filename = filename,
            SOLVEROptions...
        )
    else
        return _power_spectrum(M, b, op, ω_list; 
            solver = solver, 
            verbose = verbose,
            filename = filename,
            SOLVEROptions...
        )
    end
end

@noinline function _power_spectrum(
        M::AbstractHEOMMatrix, 
        b::AbstractVector, 
        op,
        ω_list::AbstractVector; 
        solver=UMFPACKFactorization(), 
        verbose::Bool = true,
        filename::String = "",
        SOLVEROptions...
    )

    SAVE::Bool = (filename != "")

    Size, = size(M)
    I_total = sparse(I, Size, Size)
    I_heom  = sparse(I, M.N, M.N)

    # equal to : transpose(sparse(vec(system_identity_matrix)))
    I_dual_vec = transpose(sparsevec([1 + n * (M.dim + 1) for n in 0:(M.dim - 1)], ones(M.dim), M.sup_dim))

    # operator for calculating two-time correlation functions in frequency domain
    A_normal = kron(I_heom, spre(op))
    A_dagger = kron(I_heom, spre(op'))
    local B::Vector{ComplexF64} = -1 * A_normal * b

    Length = length(ω_list)
    Sω = Vector{Float64}(undef, Length)

    if verbose
        print("Calculating spectrum for bosonic systems...\n")
        flush(stdout)
        prog = Progress(Length; start=1, desc="Progress : ", PROGBAR_OPTIONS...)
    end
    Iω   = 1im * ω_list[1] * I_total
    prob = init(LinearProblem(M.data - Iω, B), solver, SOLVEROptions...)
    sol = solve(prob)
    @inbounds for (i, ω) in enumerate(ω_list)
        if i > 1            
            Iω  = 1im * ω * I_total
            try 
                sol = solve(set_A(sol.cache, M.data - Iω))
            catch e
                if isa(e, ArgumentError)
                    prob = init(LinearProblem(M.data - Iω, B), solver, SOLVEROptions...)
                    sol = solve(prob)
                else
                    throw(e)
                end
            end
        end

        # trace over the Hilbert space of system (expectation value)
        Sω[i] = real(I_dual_vec * (A_dagger * sol.u)[1:(M.sup_dim)])

        if SAVE
            open(filename, "a") do file
                write(file, "$(psd[i]),\n")
            end
        end

        if verbose
            next!(prog)
        end 
    end
    if verbose
        println("[DONE]")
    end

    return Sω
end

@noinline function _density_of_states(
        M::AbstractHEOMMatrix, 
        b::AbstractVector, 
        op,
        ω_list::AbstractVector; 
        solver=UMFPACKFactorization(), 
        verbose::Bool = true,
        filename::String = "",
        SOLVEROptions...
    )

    SAVE::Bool = (filename != "")

    Size, = size(M)
    I_total = sparse(I, Size, Size)
    I_heom  = sparse(I, M.N, M.N)

    # equal to : transpose(sparse(vec(system_identity_matrix)))
    I_dual_vec = transpose(sparsevec([1 + n * (M.dim + 1) for n in 0:(M.dim - 1)], ones(M.dim), M.sup_dim))

    # operators for calculating two-time correlation functions in frequency domain
    A_normal = kron(I_heom, spre(op))
    A_dagger = kron(I_heom, spre(op'))
    local B_a ::Vector{ComplexF64} = -1 * A_normal * b
    local B_ad::Vector{ComplexF64} = -1 * A_dagger * b

    Length = length(ω_list)
    Aω = Vector{Float64}(undef, Length)

    if verbose
        print("Calculating spectrum for fermionic systems...\n")
        flush(stdout)
        prog = Progress(Length; start=1, desc="Progress : ", PROGBAR_OPTIONS...)
    end
    Iω = 1im * ω_list[1] * I_total
    prob_a  = init(LinearProblem(M.data - Iω, B_a),  solver, SOLVEROptions...)
    prob_ad = init(LinearProblem(M.data + Iω, B_ad), solver, SOLVEROptions...)
    sol_a  = solve(prob_a)
    sol_ad = solve(prob_ad)
    @inbounds for (i, ω) in enumerate(ω_list)
        if i > 1
            Iω = 1im * ω * I_total
            try 
                sol_a = solve(set_A(sol_a.cache, M.data - Iω))
            catch e
                if isa(e, ArgumentError)
                    prob_a = init(LinearProblem(M.data - Iω, B_a), solver, SOLVEROptions...)
                    sol_a  = solve(prob_a)
                else
                    throw(e)
                end
            end
            try
                sol_ad = solve(set_A(sol_ad.cache, M.data + Iω))
            catch e
                if isa(e, ArgumentError)
                    prob_ad = init(LinearProblem(M.data + Iω, B_ad), solver, SOLVEROptions...)
                    sol_ad  = solve(prob_ad)
                else
                    throw(e)
                end
            end
        end
        Cω_ad_a = A_dagger * sol_a.u
        Cω_a_ad = A_normal * sol_ad.u
        
        # trace over the Hilbert space of system (expectation value)
        Aω[i] = real(I_dual_vec * Cω_ad_a[1:(M.sup_dim)]) + real(I_dual_vec * Cω_a_ad[1:(M.sup_dim)])

        if SAVE
            open(filename, "a") do file
                write(file, "$(dos[i]),\n")
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
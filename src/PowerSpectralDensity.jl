"""
    PSD(M, ρ, ω_list, OP; [solver, progressBar, filename, SOLVEROptions...])
Calculate power spectral density.

# Parameters
- `M::AbstractHEOMMatrix` : the matrix given from HEOM model (the parity must be `:none` or `:even`.)
- `ρ::Union{AbstractMatrix, ADOs}` :  the system density matrix or the auxiliary density operators.
- `ω_list::AbstractVector` : the specific frequency points to solve.
- `OP::AbstractMatrix` : The system operator for the two-time correlation function in frequency domain.
- `solver` : solver in package `LinearSolve.jl`. Default to `UMFPACKFactorization()`.
- `progressBar::Bool` : Display progress bar during the process or not. Defaults to `true`.
- `filename::String` : If filename was specified, the value of psd for each ω will be saved into the file during the solving process.
- `SOLVEROptions` : extra options for solver 

# Returns
- `psd::AbstractVector` : power spectral density
"""
function PSD(
        M::AbstractHEOMMatrix, 
        ρ::T, 
        ω_list::AbstractVector, 
        OP::AbstractMatrix; 
        solver=UMFPACKFactorization(), 
        progressBar::Bool = true,
        filename::String = "",
        SOLVEROptions...
    ) where T <: Union{AbstractMatrix, ADOs}

    SAVE::Bool = (filename != "")
    if SAVE && isfile(filename)
        error("FILE: $(filename) already exist.")
    end

    # check number of bosonic states
    if (M.Nb == 0)
        error("The number of bosonic states must be greater than zero, i.e., \"M.Nb > 0\".")

    # if the bath encludes fermion states, check parity
    elseif (M.parity == :odd)
        error("The parity of M must be \":none\" (bosonic) or \":even\" (mixed) bath.")
    end

    Size, = size(M)
    I_total = sparse(I, Size, Size)
    I_heom  = sparse(I, M.N, M.N)

    local b::AbstractVector

    # check ρ
    if T <: AbstractMatrix
        if size(ρ) == (M.dim, M.dim) 
            v = sparsevec(ρ)
            b = sparsevec(v.nzind, v.nzval, Size)
        else
            error("The dimension of ρ should be equal to \"($(M.dim), $(M.dim))\".")
        end

    else # ρ::ADOs
        if (M.dim != ρ.dim)
            error("The system dimension between M and ρ are not consistent.")
        end
    
        if (M.Nb != ρ.Nb)
            error("The number of bosonic states between M and ρ are not consistent.")
        end
    
        if (M.Nf != ρ.Nf)
            error("The number of fermionic states between M and ρ are not consistent.")
        end

        b = ρ.data
    end

    # check dimension of OP
    if size(OP) != (M.dim, M.dim)
        error("The dimension of OP should be equal to \"($(M.dim), $(M.dim))\".")
    end

    # equal to : transpose(sparse(vec(system_identity_matrix)))
    I_dual_vec = transpose(sparsevec([1 + n * (M.dim + 1) for n in 0:(M.dim - 1)], ones(M.dim), M.sup_dim))

    # operator for calculating two-time correlation functions in frequency domain
    C_normal = kron(I_heom, spre(OP))
    C_dagger = kron(I_heom, spre(OP'))
    local Cb::Vector{ComplexF64} = -1 * C_normal * b

    Length = length(ω_list)
    psd    = Vector{Float64}(undef, Length)
    print("Start calculating power spectral density...")
    flush(stdout)
    if progressBar
        print("\n")
        prog = Progress(Length; start=1, desc="Progress : ", PROGBAR_OPTIONS...)
    end
    @inbounds for (i, ω) in enumerate(ω_list)
        if ω == 0
            sol = solve(LinearProblem(M.data, Cb), solver, SOLVEROptions...)
        else
            Iω = 1im * ω * I_total
            sol = solve(LinearProblem(M.data - Iω, Cb), solver, SOLVEROptions...)
        end
        Cω = C_dagger * sol.u

        # trace over the Hilbert space of system (expectation value)
        psd[i] = real(I_dual_vec * Cω[1:(M.sup_dim)])

        if SAVE
            open(filename, "a") do file
                write(file, "$(psd[i]),\n")
            end
        end

        if progressBar
            next!(prog)
        end 
    end
    println("[DONE]")

    return dos
end
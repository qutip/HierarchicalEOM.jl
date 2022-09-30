# func. for solving evolution ODE
function integrate!(dρ, ρ, L, t)
    @inbounds dρ .= L * ρ
end

"""
    evolution(M, ρ0, tlist; solver, reltol, abstol, maxiters, progressBar, SOLVEROptions...)
Solve the evolution (ODE problem) using HEOM model.

# Parameters
- `M::AbstractHEOMMatrix` : the matrix given from HEOM model
- `ρ0::AbstractMatrix` : system initial state (density matrix)
- `tlist::AbstractVector` : Denote the specific time points to save the solution at, during the solving process.
- `solver` : solver in package `DifferentialEquations.jl`. Default to `DP5()`.
- `reltol::Real` : Relative tolerance in adaptive timestepping. Default to `1.0e-6`.
- `abstol::Real` : Absolute tolerance in adaptive timestepping. Default to `1.0e-8`.
- `maxiters::Real` : Maximum number of iterations before stopping. Default to `1e5`.
- `progressBar::Bool` : Display progress bar during the process or not. Defaults to `true`.
- `SOLVEROptions` : extra options for solver 

# Returns
- `ADOs_list` : The auxiliary density operators in each time point.
"""
function evolution(
        M::AbstractHEOMMatrix, 
        ρ0::AbstractMatrix, 
        tlist::AbstractVector;
        solver = DP5(),
        reltol::Real = 1.0e-6,
        abstol::Real = 1.0e-8,
        maxiters::Real = 1e5,
        progressBar::Bool = true,
        SOLVEROptions...
    )

    if size(ρ0) == (M.dim, M.dim) 
        error("The dimension of ρ0 should be equal to \"($(M.dim), $(M.dim))\".")
    end

    # setup ρ_he 
    ADOs_list::Vector{ADOs} = []

    # vectorize initial state
    ρ0 = sparse(sparsevec(ρ0))
    ρ_he::SparseVector{ComplexF64, Int64} = sparsevec(ρ0.nzind, ρ0.nzval, M.N * M.sup_dim)
    push!(ADOs_list, ADOs(ρ_he, M.dim, M.Nb, M.Nf))
    
    # setup integrator
    dt_list = diff(tlist)
    integrator = init(
        ODEProblem(integrate!, ρ_he, (tlist[1], tlist[end]), M.data),
        solver;
        reltol = reltol,
        abstol = abstol,
        maxiters = maxiters,
        SOLVEROptions...
    )
    
    # start solving ode
    print("Start solving evolution...")
    if progressBar
        print("\n")
        prog = Progress(length(tlist); start=1, desc="Progress : ", PROGBAR_OPTIONS...)
    end
    flush(stdout)
    for dt in dt_list
        step!(integrator, dt, true)
        
        # save the ADOs
        push!(ADOs_list, ADOs(integrator.u, M.dim, M.Nb, M.Nf))
    
        if progressBar
            next!(prog)
        end
    end

    println("[DONE]\n")
    flush(stdout)

    return ADOs_list
end

"""
    evolution(M, ados, tlist; solver, reltol, abstol, maxiters, progressBar, SOLVEROptions...)
Solve the evolution (ODE problem) using HEOM model.

# Parameters
- `M::AbstractHEOMMatrix` : the matrix given from HEOM model
- `ados::ADOs` : initial auxiliary density operators
- `tlist::AbstractVector` : Denote the specific time points to save the solution at, during the solving process.
- `solver` : solver in package `DifferentialEquations.jl`. Default to `DP5()`.
- `reltol::Real` : Relative tolerance in adaptive timestepping. Default to `1.0e-6`.
- `abstol::Real` : Absolute tolerance in adaptive timestepping. Default to `1.0e-8`.
- `maxiters::Real` : Maximum number of iterations before stopping. Default to `1e5`.
- `progressBar::Bool` : Display progress bar during the process or not. Defaults to `true`.
- `SOLVEROptions` : extra options for solver 

# Returns
- `ADOs_list` : The auxiliary density operators in each time point.
"""
function evolution(
        M::AbstractHEOMMatrix, 
        ados::ADOs, 
        tlist::AbstractVector;
        solver = DP5(),
        reltol::Real = 1.0e-6,
        abstol::Real = 1.0e-8,
        maxiters::Real = 1e5,
        progressBar::Bool = true,
        SOLVEROptions...
    )

    if (M.dim != ados.dim)
        error("The system dimension between M and ados are not consistent.")
    end

    if (M.Nb != ados.Nb)
        error("The number of bosonic states between M and ados are not consistent.")
    end

    if (M.Nf != ados.Nf)
        error("The number of fermionic states between M and ados are not consistent.")
    end

    ADOs_list::Vector{ADOs} = [ados]
    
    # setup integrator
    dt_list = diff(tlist)
    integrator = init(
        ODEProblem(integrate!, ados.data, (tlist[1], tlist[end]), M.data),
        solver;
        reltol = reltol,
        abstol = abstol,
        maxiters = maxiters,
        SOLVEROptions...
    )
    
    # start solving ode
    print("Start solving evolution...")
    if progressBar
        print("\n")
        prog = Progress(length(tlist); start=1, desc="Progress : ", PROGBAR_OPTIONS...)
    end
    flush(stdout)
    for dt in dt_list
        step!(integrator, dt, true)
        
        # save the ADOs
        push!(ADOs_list, ADOs(integrator.u, M.dim, M.Nb, M.Nf))
    
        if progressBar
            next!(prog)
        end
    end

    println("[DONE]\n")
    flush(stdout)

    return ADOs_list
end
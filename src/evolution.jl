# func. for solving evolution ODE
function _hierarchy!(dρ, ρ, L, t)
    @inbounds dρ .= L * ρ
end

"""
    evolution(M, ρ0, tlist; solver, reltol, abstol, maxiters, save_everystep, progressBar, SOLVEROptions...)
Solve the time evolution for auxiliary density operators.

# Parameters
- `M::AbstractHEOMMatrix` : the matrix given from HEOM model
- `ρ0::AbstractMatrix` : system initial state (density matrix)
- `tlist::AbstractVector` : Denote the specific time points to save the solution at, during the solving process.
- `solver` : solver in package `DifferentialEquations.jl`. Default to `DP5()`.
- `reltol::Real` : Relative tolerance in adaptive timestepping. Default to `1.0e-6`.
- `abstol::Real` : Absolute tolerance in adaptive timestepping. Default to `1.0e-8`.
- `maxiters::Real` : Maximum number of iterations before stopping. Default to `1e5`.
- `save_everystep::Bool` : Saves the result at every step. Defaults to `false`.
- `progressBar::Bool` : Display progress bar during the process or not. Defaults to `true`.
- `SOLVEROptions` : extra options for solver

For more details about solvers and extra options, please refer to [`DifferentialEquations.jl`](https://diffeq.sciml.ai/stable/)

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
        save_everystep::Bool=false,
        progressBar::Bool = true,
        SOLVEROptions...
    )

    if size(ρ0) != (M.dim, M.dim) 
        error("The dimension of ρ0 should be equal to \"($(M.dim), $(M.dim))\".")
    end

    # vectorize initial state
    ρ0   = sparse(sparsevec(ρ0))
    ρ_he = sparsevec(ρ0.nzind, ρ0.nzval, M.N * M.sup_dim)
    ados = ADOs(ρ_he, M.Nb, M.Nf)
    
    return evolution(M, ados, tlist;
        solver = solver,
        reltol = reltol,
        abstol = abstol,
        maxiters = maxiters,
        save_everystep = save_everystep,
        progressBar = progressBar,
        SOLVEROptions...
    )
end

"""
    evolution(M, ados, tlist; solver, reltol, abstol, maxiters, save_everystep, progressBar, SOLVEROptions...)
Solve the time evolution for auxiliary density operators.

# Parameters
- `M::AbstractHEOMMatrix` : the matrix given from HEOM model
- `ados::ADOs` : initial auxiliary density operators
- `tlist::AbstractVector` : Denote the specific time points to save the solution at, during the solving process.
- `solver` : solver in package `DifferentialEquations.jl`. Default to `DP5()`.
- `reltol::Real` : Relative tolerance in adaptive timestepping. Default to `1.0e-6`.
- `abstol::Real` : Absolute tolerance in adaptive timestepping. Default to `1.0e-8`.
- `maxiters::Real` : Maximum number of iterations before stopping. Default to `1e5`.
- `save_everystep::Bool` : Saves the result at every step. Defaults to `false`.
- `progressBar::Bool` : Display progress bar during the process or not. Defaults to `true`.
- `SOLVEROptions` : extra options for solver

For more details about solvers and extra options, please refer to [`DifferentialEquations.jl`](https://diffeq.sciml.ai/stable/)

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
        save_everystep::Bool=false,
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
        ODEProblem(_hierarchy!, ados.data, (tlist[1], tlist[end]), M.data),
        solver;
        reltol = reltol,
        abstol = abstol,
        maxiters = maxiters,
        save_everystep = save_everystep,
        SOLVEROptions...
    )
    
    # start solving ode
    print("Solving time evolution for auxiliary density operators...")
    if progressBar
        print("\n")
        prog = Progress(length(tlist); start=1, desc="Progress : ", PROGBAR_OPTIONS...)
    end
    flush(stdout)
    for dt in dt_list
        step!(integrator, dt, true)
        
        # save the ADOs
        push!(ADOs_list, ADOs(copy(integrator.u), M.Nb, M.Nf))
    
        if progressBar
            next!(prog)
        end
    end

    GC.gc()  # clean the garbage collector
    println("[DONE]\n")
    flush(stdout)

    return ADOs_list
end
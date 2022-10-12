"""
    evolution(M, ρ0, tlist; solver, reltol, abstol, maxiters, save_everystep, verbose, filename, SOLVEROptions...)
Solve the time evolution for auxiliary density operators with initial state is given in the type of density-matrix (`ρ0`).

# Parameters
- `M::AbstractHEOMMatrix` : the matrix given from HEOM model
- `ρ0` : system initial state (density matrix)
- `tlist::AbstractVector` : Denote the specific time points to save the solution at, during the solving process.
- `solver` : solver in package `DifferentialEquations.jl`. Default to `DP5()`.
- `reltol::Real` : Relative tolerance in adaptive timestepping. Default to `1.0e-6`.
- `abstol::Real` : Absolute tolerance in adaptive timestepping. Default to `1.0e-8`.
- `maxiters::Real` : Maximum number of iterations before stopping. Default to `1e5`.
- `save_everystep::Bool` : Saves the result at every step. Defaults to `false`.
- `verbose::Bool` : To display verbose output and progress bar during the process or not. Defaults to `true`.
- `filename::String` : If filename was specified, the ADOs at each time point will be saved into the JLD2 file during the solving process.
- `SOLVEROptions` : extra options for solver

For more details about solvers and extra options, please refer to [`DifferentialEquations.jl`](https://diffeq.sciml.ai/stable/)

# Returns
- `ADOs_list` : The auxiliary density operators in each time point.
"""
function evolution(
        M::AbstractHEOMMatrix, 
        ρ0, 
        tlist::AbstractVector;
        solver = DP5(),
        reltol::Real = 1.0e-6,
        abstol::Real = 1.0e-8,
        maxiters::Real = 1e5,
        save_everystep::Bool=false,
        verbose::Bool = true,
        filename::String = "",
        SOLVEROptions...
    )

    if !isValidMatrixType(ρ0, M.dim)
        error("Invalid matrix \"ρ0\".")
    end

    # vectorize initial state
    ρ1   = sparse(sparsevec(ρ0))
    ados = ADOs(
        sparsevec(ρ1.nzind, ρ1.nzval, M.N * M.sup_dim), 
        M.Nb, 
        M.Nf
    )
    
    return evolution(M, ados, tlist;
        solver = solver,
        reltol = reltol,
        abstol = abstol,
        maxiters = maxiters,
        save_everystep = save_everystep,
        verbose = verbose,
        filename = filename,
        SOLVEROptions...
    )
end

"""
    evolution(M, ados, tlist; solver, reltol, abstol, maxiters, save_everystep, verbose, filename, SOLVEROptions...)
Solve the time evolution for auxiliary density operators with initial state is given in the type of `ADOs`.

# Parameters
- `M::AbstractHEOMMatrix` : the matrix given from HEOM model
- `ados::ADOs` : initial auxiliary density operators
- `tlist::AbstractVector` : Denote the specific time points to save the solution at, during the solving process.
- `solver` : solver in package `DifferentialEquations.jl`. Default to `DP5()`.
- `reltol::Real` : Relative tolerance in adaptive timestepping. Default to `1.0e-6`.
- `abstol::Real` : Absolute tolerance in adaptive timestepping. Default to `1.0e-8`.
- `maxiters::Real` : Maximum number of iterations before stopping. Default to `1e5`.
- `save_everystep::Bool` : Saves the result at every step. Defaults to `false`.
- `verbose::Bool` : To display verbose output and progress bar during the process or not. Defaults to `true`.
- `filename::String` : If filename was specified, the ADOs at each time point will be saved into the JLD2 file during the solving process.
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
        verbose::Bool = true,
        filename::String = "",
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

    SAVE::Bool = (filename != "")
    if SAVE && isfile(filename)
        error("FILE: $(filename) already exist.")
    end

    ADOs_list::Vector{ADOs} = [ados]
    if SAVE
        jldopen(filename, "a") do file
            file[string(tlist[1])] = ados
        end
    end
    
    # setup ode function
    S, = size(M)
    hierarchy = ODEFunction(_hierarchy!; jac_prototype = spzeros(ComplexF64, S, S))

    # setup integrator
    integrator = init(
        ODEProblem(hierarchy, Vector(ados.data), (tlist[1], tlist[end]), M.data),
        solver;
        reltol = reltol,
        abstol = abstol,
        maxiters = maxiters,
        save_everystep = save_everystep,
        SOLVEROptions...
    )
    
    # start solving ode
    if verbose
        print("Solving time evolution for auxiliary density operators...\n")
        flush(stdout)
        prog = Progress(length(tlist); start=1, desc="Progress : ", PROGBAR_OPTIONS...)
    end
    idx = 1
    dt_list = diff(tlist)
    for dt in dt_list
        idx += 1
        step!(integrator, dt, true)
        
        # save the ADOs
        ados = ADOs(copy(integrator.u), M.Nb, M.Nf)
        push!(ADOs_list, ados)
        if SAVE
            jldopen(filename, "a") do file
                file[string(tlist[idx])] = ados
            end
        end
    
        if verbose
            next!(prog)
        end
    end

    GC.gc()  # clean the garbage collector
    if verbose
        println("[DONE]\n")
        flush(stdout)
    end

    return ADOs_list
end
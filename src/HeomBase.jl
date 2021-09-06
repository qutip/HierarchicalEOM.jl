abstract type AbstractHEOMMatrix end
size(A::AbstractHEOMMatrix) = size(A.data)

include("ados.jl")

spre(q::AbstractMatrix)        = kron(Matrix(I, size(q)[1], size(q)[1]), q)
spre(q::AbstractOperator)      = sparse(kron(sparse(I, size(q)[1], size(q)[1]), q.data))
spre(q::AbstractSparseMatrix)  = sparse(kron(sparse(I, size(q)[1], size(q)[1]), q))
spost(q::AbstractMatrix)       = kron(transpose(q), Matrix(I, size(q)[1], size(q)[1]))
spost(q::AbstractOperator)     = sparse(kron(transpose(q.data), sparse(I, size(q)[1], size(q)[1])))
spost(q::AbstractSparseMatrix) = sparse(kron(transpose(q), sparse(I, size(q)[1], size(q)[1])))

# generate liouvillian matrix
function liouvillian(Hsys, Jump_Ops::Vector=[], progressBar::Bool=true)
        
    N, = size(Hsys)

    L = -1im * (spre(Hsys) - spost(Hsys))
    if progressBar
        prog = Progress(length(Jump_Ops) + 1, start=1; desc="Construct Liouvillian     : ", PROGBAR_OPTIONS...)
    end
    for J in Jump_Ops
        L += spre(J) * spost(J') - 0.5 * (spre(J' * J) + spost(J' * J))
        if progressBar
            next!(prog)
        end
    end
    return L
end

# func. for solving evolution ODE
function hierachy!(dρ, ρ, L, t)
    dρ .= L * ρ
end

"""
# `evolution(M, ρ0, tlist; [solver, returnAdos, progressBar, SOLVEROptions...])`
Solve the evolution (ODE problem) using HEOM model.

## Parameters
- `M::AbstractHEOMMatrix` : the matrix given from HEOM model
- `ρ0::Union{AbstractMatrix, AbstractOperator}` : initial state (density matrix)
- `tlist::AbstractVector` : Denote the specific time points to save the solution at, during the solving process.
- `solver` : solver in package `DifferentialEquations.jl`. Default to `Tsit5()`.
- `returnAdos::Bool` : Decide to return the entire Ados vector or not. Default as `false`.
- `progressBar::Bool` : Display progress bar during the process or not. Defaults to `true`.
- `SOLVEROptions` : extra options for solver 

## Returns
- `ρ_list` : The reduced density matrices in each time point.
- `Ados_vec` : The Ados vector in each time point (only return when `returnAdos=true`).
"""
function evolution(
        M::AbstractHEOMMatrix, 
        ρ0::Union{AbstractMatrix, AbstractOperator}, 
        tlist::AbstractVector;
        solver=Tsit5(),
        returnAdos::Bool=false,
        progressBar::Bool=true,
        SOLVEROptions...
    )

    (N1, N2) = size(ρ0)
    if (N1 != M.N_sys) || (N2 != M.N_sys)
        error("The size of initial state ρ0 is incorrect.")
    end

    # setup ρ_he and ρlist
    ρlist::Vector{SparseMatrixCSC{ComplexF64, Int64}} = []
    ρ_he::SparseVector{ComplexF64, Int64} = spzeros(M.N_he * M.sup_dim)
    if typeof(ρ0) <: AbstractMatrix
        push!(ρlist, ρ0)
        ρ_he[1:(M.sup_dim)] .= vec(ρ0)
    else
        push!(ρlist, ρ0.data)
        ρ_he[1:(M.sup_dim)] .= vec(ρ0.data)
    end

    if returnAdos
        Ados = [ρ_he]
    end
    
    # setup integrator
    dt_list = diff(tlist)
    integrator = init(
            ODEProblem(hierachy!, ρ_he, (tlist[1], tlist[end]), M.data),
            solver;
            SOLVEROptions...
    )
    
    # start solving ode
    print("Start solving hierachy equation of motions...")
    if progressBar
        print("\n")
        prog = Progress(length(tlist); start=1, desc="Progress : ", PROGBAR_OPTIONS...)
    end
    flush(stdout)
    for dt in dt_list
        step!(integrator, dt, true)
        
        # save the reduced density matrix and current ados vector (if \'returnAdos\' is specified)
        push!(ρlist, reshape(integrator.u[1:(M.sup_dim)], M.N_sys, M.N_sys))
        if returnAdos
            push!(Ados, integrator.u)
        end
    
        if progressBar
            next!(prog)
        end
    end

    println("[DONE]\n")
    flush(stdout)

    if returnAdos
        return ρlist, Ados
    else
        return ρlist
    end
end
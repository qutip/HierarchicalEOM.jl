export M_Boson

@doc raw"""
    struct M_Boson <: AbstractHEOMLSMatrix
HEOM Liouvillian superoperator matrix for bosonic bath

# Fields
- `data<:AbstractSciMLOperator` : the matrix of HEOM Liouvillian superoperator
- `tier` : the tier (cutoff level) for the bosonic hierarchy
- `dimensions` : the dimension list of the coupling operator (should be equal to the system dimensions).
- `N` : the number of total ADOs
- `sup_dim` : the dimension of system superoperator
- `parity` : the parity label of the operator which HEOMLS is acting on (usually `EVEN`, only set as `ODD` for calculating spectrum of fermionic system).
- `bath::Vector{BosonBath}` : the vector which stores all `BosonBath` objects
- `hierarchy::HierarchyDict`: the object which contains all dictionaries for boson-bath-ADOs hierarchy.

!!! note "`dims` property"
    For a given `M::M_Boson`, `M.dims` or `getproperty(M, :dims)` returns its `dimensions` in the type of integer-vector.
"""
struct M_Boson{T<:AbstractSciMLOperator} <: AbstractHEOMLSMatrix{T}
    data::T
    tier::Int
    dimensions::Dimensions
    N::Int
    sup_dim::Int
    parity::AbstractParity
    bath::Vector{BosonBath}
    hierarchy::HierarchyDict
end

function M_Boson(
    Hsys::QuantumObject,
    tier::Int,
    Bath::BosonBath,
    parity::AbstractParity = EVEN;
    threshold::Real = 0.0,
    verbose::Bool = true,
)
    return M_Boson(Hsys, tier, [Bath], parity, threshold = threshold, verbose = verbose)
end

@doc raw"""
    M_Boson(Hsys, tier, Bath, parity=EVEN; threshold=0.0, verbose=true)
Generate the boson-type HEOM Liouvillian superoperator matrix

# Parameters
- `Hsys` : The time-independent system Hamiltonian or Liouvillian
- `tier::Int` : the tier (cutoff level) for the bosonic bath
- `Bath::Vector{BosonBath}` : objects for different bosonic baths
- `parity::AbstractParity` : the parity label of the operator which HEOMLS is acting on (usually `EVEN`, only set as `ODD` for calculating spectrum of fermionic system).
- `threshold::Real` : The threshold of the importance value (see Ref. [1]). Defaults to `0.0`.
- `verbose::Bool` : To display verbose output and progress bar during the process or not. Defaults to `true`.

Note that the parity only need to be set as `ODD` when the system contains fermionic systems and you need to calculate the spectrum (density of states) of it.

[1] [Phys. Rev. B 88, 235426 (2013)](https://doi.org/10.1103/PhysRevB.88.235426)
"""
@noinline function M_Boson(
    Hsys::QuantumObject,
    tier::Int,
    Bath::Vector{BosonBath},
    parity::AbstractParity = EVEN;
    threshold::Real = 0.0,
    verbose::Bool = true,
)

    # check for system dimension
    _Hsys = HandleMatrixType(Hsys, "Hsys (system Hamiltonian or Liouvillian)")
    sup_dim = prod(_Hsys.dimensions)^2
    I_sup = sparse(one(ComplexF64) * I, sup_dim, sup_dim)

    # the Liouvillian operator for free Hamiltonian term
    Lsys = minus_i_L_op(_Hsys)

    # bosonic bath
    if verbose && (threshold > 0.0)
        print("Checking the importance value for each ADOs...")
        flush(stdout)
    end
    Nado, baths, hierarchy = genBathHierarchy(Bath, tier, _Hsys.dimensions, threshold = threshold)
    idx2nvec = hierarchy.idx2nvec
    nvec2idx = hierarchy.nvec2idx
    if verbose && (threshold > 0.0)
        println("[DONE]")
        flush(stdout)
    end

    # time dependent scalar operators
    td_scalars = ScalarOperator[]
    FunctionFieldType = Union{bosonInputFunction,bosonOutputFunctionLeft,bosonOutputFunctionRight}
    for b in baths
        (b isa FunctionFieldType) && append!(td_scalars, b.η)
    end

    # start to construct the matrix
    Nthread = nthreads()
    λ0 = ScalarOperator(0) # this is just a key for conventional (time-independent) HEOMLS
    L_row = Dict(λ0 => Vector{Int}[Int[] for _ in 1:Nthread])
    L_col = Dict(λ0 => Vector{Int}[Int[] for _ in 1:Nthread])
    L_val = Dict(λ0 => Vector{ComplexF64}[ComplexF64[] for _ in 1:Nthread])
    for λ in td_scalars
        L_row[λ] = Vector{Int}[Int[] for _ in 1:Nthread]
        L_col[λ] = Vector{Int}[Int[] for _ in 1:Nthread]
        L_val[λ] = Vector{ComplexF64}[ComplexF64[] for _ in 1:Nthread]
    end

    if verbose
        println("Preparing block matrices for HEOM Liouvillian superoperator (using $(Nthread) threads)...")
        flush(stdout)
        prog = ProgressBar(Nado)
    end
    @threads for idx in 1:Nado
        tID = threadid()

        # boson (current level) superoperator
        nvec = idx2nvec[idx]
        if nvec.level >= 1
            sum_γ = bath_sum_γ(nvec, baths)
            op = Lsys - sum_γ * I_sup
        else
            op = Lsys
        end
        add_operator!(op, L_row[λ0][tID], L_col[λ0][tID], L_val[λ0][tID], Nado, idx, idx)

        # connect to bosonic (n+1)th- & (n-1)th- level superoperator
        mode = 0
        nvec_neigh = copy(nvec)
        for bB in baths
            for k in 1:bB.Nterm
                mode += 1
                n_k = nvec[mode]

                # connect to bosonic (n-1)th-level superoperator
                if n_k > 0
                    Nvec_minus!(nvec_neigh, mode)
                    if haskey(nvec2idx, nvec_neigh)
                        idx_neigh = nvec2idx[nvec_neigh]
                        op = minus_i_D_op(bB, k, n_k)
                        if bB isa FunctionFieldType
                            λ = bB.η[k]
                        else
                            λ = λ0
                        end
                        add_operator!(op, L_row[λ][tID], L_col[λ][tID], L_val[λ][tID], Nado, idx, idx_neigh)
                    end
                    Nvec_plus!(nvec_neigh, mode)
                end

                # connect to bosonic (n+1)th-level superoperator
                if nvec.level < tier
                    Nvec_plus!(nvec_neigh, mode)
                    if haskey(nvec2idx, nvec_neigh)
                        idx_neigh = nvec2idx[nvec_neigh]
                        op = minus_i_B_op(bB)
                        add_operator!(op, L_row[λ0][tID], L_col[λ0][tID], L_val[λ0][tID], Nado, idx, idx_neigh)
                    end
                    Nvec_minus!(nvec_neigh, mode)
                end
            end
        end
        if verbose
            next!(prog) # trigger a progress bar update
        end
    end

    if verbose
        print("Constructing matrix...")
        flush(stdout)
    end

    # conventional HEOMLS (time independent)
    L_he = MatrixOperator(
        sparse(
            reduce(vcat, L_row[λ0]),
            reduce(vcat, L_col[λ0]),
            reduce(vcat, L_val[λ0]),
            Nado * sup_dim,
            Nado * sup_dim,
        ),
    )

    # Input/Output HEOMLS (time dependent)
    if !isempty(td_scalars)
        L_he += mapreduce(
            λ ->
                λ * MatrixOperator(
                    sparse(
                        reduce(vcat, L_row[λ]),
                        reduce(vcat, L_col[λ]),
                        reduce(vcat, L_val[λ]),
                        Nado * sup_dim,
                        Nado * sup_dim,
                    ),
                ),
            +,
            td_scalars,
        )
    end

    if verbose
        println("[DONE]")
        flush(stdout)
    end
    return M_Boson(L_he, tier, _Hsys.dimensions, Nado, sup_dim, parity, Bath, hierarchy)
end

_getBtier(M::M_Boson) = M.tier
_getFtier(M::M_Boson) = 0

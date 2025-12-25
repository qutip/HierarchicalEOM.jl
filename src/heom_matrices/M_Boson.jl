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
    assemble::Union{Val,Symbol} = Val(:full),
    verbose::Bool = true,
)
    return M_Boson(Hsys, tier, [Bath], parity; threshold, assemble, verbose)
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
- `assemble::Union{Val,Symbol}` : Methods to assemble (evaluate lazy operations) the HEOMLS. Accepted values are `:full`, `:combine`, or `:none`. Defaults to `Val(:full)`.
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
    assemble::Union{Val,Symbol} = Val(:full),
    verbose::Bool = true,
)
    _Hsys = HandleMatrixType(Hsys, "Hsys (system Hamiltonian or Liouvillian)") # Checking input type first

    assemble_method = makeVal(assemble)
    check_assemble_method(assemble_method)

    # bosonic bath
    if verbose && (threshold > 0.0)
        print("Checking the importance value for each ADOs...")
        flush(stdout)
    end
    Nado, baths, hierarchy = genBathHierarchy(Bath, tier, _Hsys.dimensions; threshold)
    idx2nvec = hierarchy.idx2nvec
    nvec2idx = hierarchy.nvec2idx
    if verbose && (threshold > 0.0)
        println("[DONE]")
        flush(stdout)
    end

    minus_γ_term = zeros(ComplexF64, Nado)

    # stores position and prefix value for each Boson superoperators in HEOM Liouville space using sparse COO format
    B_terms = [HEOMSparseStructure(bB, Nado) for bB in baths]

    if verbose
        println("Preparing HEOM Liouvillian sparsity structure...")
        flush(stdout)
        progr = Progress(Nado; enabled = verbose, desc = "[M_Boson] ", QuantumToolbox.settings.ProgressMeterKWARGS...)
    end
    for idx in 1:Nado
        # boson (current level) superoperator
        nvec = idx2nvec[idx]
        if nvec.level >= 1
            minus_γ_term[idx] = -bath_sum_γ(nvec, baths)
        end

        # connect to bosonic (n+1)th- & (n-1)th- level superoperator
        mode = 0
        nvec_neigh = copy(nvec)
        for (b_term, bB) in zip(B_terms, baths)
            for k in 1:bB.Nterm
                mode += 1
                n_k = nvec[mode]

                # connect to bosonic (n-1)th-level superoperator
                if n_k > 0
                    Nvec_minus!(nvec_neigh, mode)
                    if (threshold == 0.0) || haskey(nvec2idx, nvec_neigh)
                        idx_neigh = nvec2idx[nvec_neigh]
                        minus_i_D_op!(b_term, idx, idx_neigh, bB, k, n_k)
                    end
                    Nvec_plus!(nvec_neigh, mode)
                end

                # connect to bosonic (n+1)th-level superoperator
                if nvec.level < tier
                    Nvec_plus!(nvec_neigh, mode)
                    if (threshold == 0.0) || haskey(nvec2idx, nvec_neigh)
                        idx_neigh = nvec2idx[nvec_neigh]
                        minus_i_B_op!(b_term, idx, idx_neigh, bB)
                    end
                    Nvec_minus!(nvec_neigh, mode)
                end
            end
        end
        verbose && next!(progr) # trigger a progress bar update
    end

    # Create SciML lazy HEOM Liouvillian superoperator
    sup_dim = prod(_Hsys.dimensions)^2
    L_t_indep = kron(MatrixOperator(Eye(Nado)), minus_i_L_op(_Hsys)) # the Liouvillian operator for free Hamiltonian term
    L_t_indep += kron(MatrixOperator(spdiagm(minus_γ_term)), Eye(sup_dim)) # minus sum γ terms

    # Superoperator cross level terms
    for b_term in B_terms
        for op in HEOMSparseStructureFieldNames
            b_coo = getfield(b_term, op)
            b_coo isa Nothing || (L_t_indep += _gen_HEOMLS_term(b_coo))
        end
    end
    if verbose
        println("[DONE]")
        flush(stdout)
    end

    L_heom = assemble_HEOMLS_terms(L_t_indep, assemble_method, verbose)[1]
    return M_Boson(L_heom, tier, _Hsys.dimensions, Nado, sup_dim, parity, Bath, hierarchy)
end

_getBtier(M::M_Boson) = M.tier
_getFtier(M::M_Boson) = 0

export M_Fermion

@doc raw"""
    struct M_Fermion <: AbstractHEOMLSMatrix
HEOM Liouvillian superoperator matrix for fermionic bath

# Fields
- `data<:AbstractSciMLOperator` : the matrix of HEOM Liouvillian superoperator
- `tier` : the tier (cutoff level) for the fermionic hierarchy
- `dimensions` : the dimension list of the coupling operator (should be equal to the system dimensions).
- `N` : the number of total ADOs
- `sup_dim` : the dimension of system superoperator
- `parity` : the parity label of the operator which HEOMLS is acting on (usually `EVEN`, only set as `ODD` for calculating spectrum of fermionic system).
- `bath::Vector{FermionBath}` : the vector which stores all `FermionBath` objects
- `hierarchy::HierarchyDict`: the object which contains all dictionaries for fermion-bath-ADOs hierarchy.

!!! note "`dims` property"
    For a given `M::M_Fermion`, `M.dims` or `getproperty(M, :dims)` returns its `dimensions` in the type of integer-vector.
"""
struct M_Fermion{T<:AbstractSciMLOperator} <: AbstractHEOMLSMatrix{T}
    data::T
    tier::Int
    dimensions::Dimensions
    N::Int
    sup_dim::Int
    parity::AbstractParity
    bath::Vector{FermionBath}
    hierarchy::HierarchyDict
end

function M_Fermion(
    Hsys::QuantumObject,
    tier::Int,
    Bath::FermionBath,
    parity::AbstractParity = EVEN;
    threshold::Real = 0.0,
    concretize::Union{Val,Bool} = Val(true),
    verbose::Bool = true,
)
    return M_Fermion(Hsys, tier, [Bath], parity; threshold, concretize, verbose)
end

@doc raw"""
    M_Fermion(Hsys, tier, Bath, parity=EVEN; threshold=0.0, verbose=true)
Generate the fermion-type HEOM Liouvillian superoperator matrix

# Parameters
- `Hsys` : The time-independent system Hamiltonian or Liouvillian
- `tier::Int` : the tier (cutoff level) for the fermionic bath
- `Bath::Vector{FermionBath}` : objects for different fermionic baths
- `parity::AbstractParity` : the parity label of the operator which HEOMLS is acting on (usually `EVEN`, only set as `ODD` for calculating spectrum of fermionic system).
- `threshold::Real` : The threshold of the importance value (see Ref. [1]). Defaults to `0.0`.
- `concretize::Union{Val,Bool}` : Whether to concretize the HEOMLS to a single sparse matrix. Defaults to `Val(true)`.
- `verbose::Bool` : To display verbose output and progress bar during the process or not. Defaults to `true`.

[1] [Phys. Rev. B 88, 235426 (2013)](https://doi.org/10.1103/PhysRevB.88.235426)
"""
@noinline function M_Fermion(
    Hsys::QuantumObject,
    tier::Int,
    Bath::Vector{FermionBath},
    parity::AbstractParity = EVEN;
    threshold::Real = 0.0,
    concretize::Union{Val,Bool} = Val(true),
    verbose::Bool = true,
)
    _Hsys = HandleMatrixType(Hsys, "Hsys (system Hamiltonian or Liouvillian)") # Checking input type first

    # fermionic bath
    if verbose && (threshold > 0.0)
        print("Checking the importance value for each ADOs...")
        flush(stdout)
    end
    Nado, baths, hierarchy = genBathHierarchy(Bath, tier, Hsys.dimensions; threshold)
    idx2nvec = hierarchy.idx2nvec
    nvec2idx = hierarchy.nvec2idx
    if verbose && (threshold > 0.0)
        println("[DONE]")
        flush(stdout)
    end

    γ_term = Vector{ComplexF64}(undef, Nado)
    γ_term[1] = 0.0im

    # stores position and prefix value for each Fermion superoperators in HEOM Liouville space using sparse COO format
    F_terms = [HEOMSparseStructure(fB, Nado) for fB in baths]

    if verbose
        println("Preparing HEOM Liouvillian sparsity structure...")
        flush(stdout)
        progr = Progress(Nado; enabled = verbose, desc = "[M_Fermion] ", QuantumToolbox.settings.ProgressMeterKWARGS...)
    end
    for idx in 1:Nado
        nvec = idx2nvec[idx]
        if nvec.level > 0
            γ_term[idx] = -bath_sum_γ(nvec, baths)
        end

        mode = 0
        nvec_neigh = copy(nvec)
        for (f_term, fB) in zip(F_terms, baths)
            for k in 1:fB.Nterm
                mode += 1
                n_k = nvec[mode]

                if n_k > 0
                    Nvec_minus!(nvec_neigh, mode)
                    if (threshold == 0.0) || haskey(nvec2idx, nvec_neigh)
                        idx_neigh = nvec2idx[nvec_neigh]
                        minus_i_C_op!(f_term, idx, idx_neigh, fB, k, nvec.level, sum(nvec_neigh[1:(mode-1)]), parity)
                    end
                    Nvec_plus!(nvec_neigh, mode)
                elseif nvec.level < tier
                    Nvec_plus!(nvec_neigh, mode)
                    if (threshold == 0.0) || haskey(nvec2idx, nvec_neigh)
                        idx_neigh = nvec2idx[nvec_neigh]
                        minus_i_A_op!(f_term, idx, idx_neigh, fB, nvec.level, sum(nvec_neigh[1:(mode-1)]), parity)
                    end
                    Nvec_minus!(nvec_neigh, mode)
                end
            end

        end
        verbose && next!(progr) # trigger a progress bar update
    end
    
    # Create SciML lazy HEOM Liouvillian superoperator
    sup_dim = prod(_Hsys.dimensions)^2
    L_heom = kron(MatrixOperator(Eye(Nado)), minus_i_L_op(_Hsys)) # the Liouvillian operator for free Hamiltonian term
    L_heom += kron(MatrixOperator(spdiagm(γ_term)), Eye(sup_dim)) # ADOs sum γ terms
    
    # Super operator cross level terms
    for (f_term, fB) in zip(F_terms, baths)
        for op in fieldnames(HEOMSparseStructure)
            f_coo = getfield(f_term, op)
            f_coo isa Nothing && continue
            L_heom += kron(MatrixOperator(sparse(f_coo)), getfield(fB, op))
        end
    end
    if verbose
        println("[DONE]")
        flush(stdout)
    end
    
    if getVal(makeVal(concretize))
        if verbose
            print("Concretizing matrix...")
            flush(stdout)
        end
        data = SciMLOperators.concretize(L_heom) |> MatrixOperator
        if verbose
            println("[DONE]")
            flush(stdout)
        end
    else
        data = L_heom
    end

    return M_Fermion(data, tier, _Hsys.dimensions, Nado, sup_dim, parity, Bath, hierarchy)
end

_getBtier(M::M_Fermion) = 0
_getFtier(M::M_Fermion) = M.tier

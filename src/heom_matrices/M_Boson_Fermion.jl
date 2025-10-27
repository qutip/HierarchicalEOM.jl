export M_Boson_Fermion

@doc raw"""
    struct M_Boson_Fermion <: AbstractHEOMLSMatrix
HEOM Liouvillian superoperator matrix for mixtured (bosonic and fermionic) bath 

# Fields
- `data<:AbstractSciMLOperator` : the matrix of HEOM Liouvillian superoperator
- `Btier` : the tier (cutoff level) for bosonic hierarchy
- `Ftier` : the tier (cutoff level) for fermionic hierarchy
- `dimensions` : the dimension list of the coupling operator (should be equal to the system dimensions).
- `N` : the number of total ADOs
- `sup_dim` : the dimension of system superoperator
- `parity` : the parity label of the operator which HEOMLS is acting on (usually `EVEN`, only set as `ODD` for calculating spectrum of fermionic system).
- `Bbath::Vector{BosonBath}` : the vector which stores all `BosonBath` objects
- `Fbath::Vector{FermionBath}` : the vector which stores all `FermionBath` objects
- `hierarchy::MixHierarchyDict`: the object which contains all dictionaries for mixed-bath-ADOs hierarchy.

!!! note "`dims` property"
    For a given `M::M_Boson_Fermion`, `M.dims` or `getproperty(M, :dims)` returns its `dimensions` in the type of integer-vector.
"""
struct M_Boson_Fermion{T<:AbstractSciMLOperator} <: AbstractHEOMLSMatrix{T}
    data::T
    Btier::Int
    Ftier::Int
    dimensions::Dimensions
    N::Int
    sup_dim::Int
    parity::AbstractParity
    Bbath::Vector{BosonBath}
    Fbath::Vector{FermionBath}
    hierarchy::MixHierarchyDict
end

function M_Boson_Fermion(
    Hsys::QuantumObject,
    Btier::Int,
    Ftier::Int,
    Bbath::BosonBath,
    Fbath::FermionBath,
    parity::AbstractParity = EVEN;
    threshold::Real = 0.0,
    verbose::Bool = true,
)
    return M_Boson_Fermion(Hsys, Btier, Ftier, [Bbath], [Fbath], parity, threshold = threshold, verbose = verbose)
end

function M_Boson_Fermion(
    Hsys::QuantumObject,
    Btier::Int,
    Ftier::Int,
    Bbath::Vector{BosonBath},
    Fbath::FermionBath,
    parity::AbstractParity = EVEN;
    threshold::Real = 0.0,
    verbose::Bool = true,
)
    return M_Boson_Fermion(Hsys, Btier, Ftier, Bbath, [Fbath], parity, threshold = threshold, verbose = verbose)
end

function M_Boson_Fermion(
    Hsys::QuantumObject,
    Btier::Int,
    Ftier::Int,
    Bbath::BosonBath,
    Fbath::Vector{FermionBath},
    parity::AbstractParity = EVEN;
    threshold::Real = 0.0,
    verbose::Bool = true,
)
    return M_Boson_Fermion(Hsys, Btier, Ftier, [Bbath], Fbath, parity, threshold = threshold, verbose = verbose)
end

@doc raw"""
    M_Boson_Fermion(Hsys, Btier, Ftier, Bbath, Fbath, parity=EVEN; threshold=0.0, verbose=true)
Generate the boson-fermion-type HEOM Liouvillian superoperator matrix

# Parameters
- `Hsys` : The time-independent system Hamiltonian or Liouvillian
- `Btier::Int` : the tier (cutoff level) for the bosonic bath
- `Ftier::Int` : the tier (cutoff level) for the fermionic bath
- `Bbath::Vector{BosonBath}` : objects for different bosonic baths
- `Fbath::Vector{FermionBath}` : objects for different fermionic baths
- `parity::AbstractParity` : the parity label of the operator which HEOMLS is acting on (usually `EVEN`, only set as `ODD` for calculating spectrum of fermionic system).
- `threshold::Real` : The threshold of the importance value (see Ref. [1, 2]). Defaults to `0.0`.
- `verbose::Bool` : To display verbose output and progress bar during the process or not. Defaults to `true`.

Note that the parity only need to be set as `ODD` when the system contains fermion systems and you need to calculate the spectrum of it.

[1] [Phys. Rev. B  88, 235426 (2013)](https://doi.org/10.1103/PhysRevB.88.235426)
[2] [Phys. Rev. B 103, 235413 (2021)](https://doi.org/10.1103/PhysRevB.103.235413)
"""
@noinline function M_Boson_Fermion(
    Hsys::QuantumObject,
    Btier::Int,
    Ftier::Int,
    Bbath::Vector{BosonBath},
    Fbath::Vector{FermionBath},
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

    # check for bosonic and fermionic bath
    if verbose && (threshold > 0.0)
        print("Checking the importance value for each ADOs...")
        flush(stdout)
    end
    Nado, baths_b, baths_f, hierarchy =
        genBathHierarchy(Bbath, Fbath, Btier, Ftier, _Hsys.dimensions, threshold = threshold)
    idx2nvec = hierarchy.idx2nvec
    nvec2idx = hierarchy.nvec2idx
    if verbose && (threshold > 0.0)
        println("[DONE]")
        flush(stdout)
    end

    # start to construct the matrix
    Nthread = nthreads()
    L_row = [Int[] for _ in 1:Nthread]
    L_col = [Int[] for _ in 1:Nthread]
    L_val = [ComplexF64[] for _ in 1:Nthread]
    chnl = Channel{Tuple{Vector{Int},Vector{Int},Vector{ComplexF64}}}(Nthread)
    foreach(i -> put!(chnl, (L_row[i], L_col[i], L_val[i])), 1:Nthread)
    if verbose
        println("Preparing block matrices for HEOM Liouvillian superoperator (using $(Nthread) threads)...")
        flush(stdout)
        progr = Progress(Nado, showspeed = true, enabled = verbose, desc = "[M_Boson_Fermion] ")
    end
    @threads for idx in 1:Nado

        # boson and fermion (current level) superoperator
        sum_γ = 0.0
        nvec_b, nvec_f = idx2nvec[idx]
        if nvec_b.level >= 1
            sum_γ += bath_sum_γ(nvec_b, baths_b)
        end
        if nvec_f.level >= 1
            sum_γ += bath_sum_γ(nvec_f, baths_f)
        end
        op = Lsys - sum_γ * I_sup
        L_tuple = take!(chnl)
        add_operator!(op, L_tuple[1], L_tuple[2], L_tuple[3], Nado, idx, idx)
        put!(chnl, L_tuple)

        # connect to bosonic (n+1)th- & (n-1)th- level superoperator
        mode = 0
        nvec_neigh = copy(nvec_b)
        for bB in baths_b
            for k in 1:bB.Nterm
                mode += 1
                n_k = nvec_b[mode]

                # connect to bosonic (n-1)th-level superoperator
                if n_k > 0
                    Nvec_minus!(nvec_neigh, mode)
                    if (threshold == 0.0) || haskey(nvec2idx, (nvec_neigh, nvec_f))
                        idx_neigh = nvec2idx[(nvec_neigh, nvec_f)]
                        op = minus_i_D_op(bB, k, n_k)
                        L_tuple = take!(chnl)
                        add_operator!(op, L_tuple[1], L_tuple[2], L_tuple[3], Nado, idx, idx_neigh)
                        put!(chnl, L_tuple)
                    end
                    Nvec_plus!(nvec_neigh, mode)
                end

                # connect to bosonic (n+1)th-level superoperator
                if nvec_b.level < Btier
                    Nvec_plus!(nvec_neigh, mode)
                    if (threshold == 0.0) || haskey(nvec2idx, (nvec_neigh, nvec_f))
                        idx_neigh = nvec2idx[(nvec_neigh, nvec_f)]
                        op = minus_i_B_op(bB)
                        L_tuple = take!(chnl)
                        add_operator!(op, L_tuple[1], L_tuple[2], L_tuple[3], Nado, idx, idx_neigh)
                        put!(chnl, L_tuple)
                    end
                    Nvec_minus!(nvec_neigh, mode)
                end
            end
        end

        # connect to fermionic (n+1)th- & (n-1)th- level superoperator
        mode = 0
        nvec_neigh = copy(nvec_f)
        for fB in baths_f
            for k in 1:fB.Nterm
                mode += 1
                n_k = nvec_f[mode]

                # connect to fermionic (n-1)th-level superoperator
                if n_k > 0
                    Nvec_minus!(nvec_neigh, mode)
                    if (threshold == 0.0) || haskey(nvec2idx, (nvec_b, nvec_neigh))
                        idx_neigh = nvec2idx[(nvec_b, nvec_neigh)]
                        op = minus_i_C_op(fB, k, nvec_f.level, sum(nvec_neigh[1:(mode-1)]), parity)
                        L_tuple = take!(chnl)
                        add_operator!(op, L_tuple[1], L_tuple[2], L_tuple[3], Nado, idx, idx_neigh)
                        put!(chnl, L_tuple)
                    end
                    Nvec_plus!(nvec_neigh, mode)

                    # connect to fermionic (n+1)th-level superoperator
                elseif nvec_f.level < Ftier
                    Nvec_plus!(nvec_neigh, mode)
                    if (threshold == 0.0) || haskey(nvec2idx, (nvec_b, nvec_neigh))
                        idx_neigh = nvec2idx[(nvec_b, nvec_neigh)]
                        op = minus_i_A_op(fB, nvec_f.level, sum(nvec_neigh[1:(mode-1)]), parity)
                        L_tuple = take!(chnl)
                        add_operator!(op, L_tuple[1], L_tuple[2], L_tuple[3], Nado, idx, idx_neigh)
                        put!(chnl, L_tuple)
                    end
                    Nvec_minus!(nvec_neigh, mode)
                end
            end
        end
        if verbose
            next!(progr) # trigger a progress bar update
        end
    end
    if verbose
        print("Constructing matrix...")
        flush(stdout)
    end
    L_he = MatrixOperator(
        sparse(reduce(vcat, L_row), reduce(vcat, L_col), reduce(vcat, L_val), Nado * sup_dim, Nado * sup_dim),
    )
    if verbose
        println("[DONE]")
        flush(stdout)
    end
    return M_Boson_Fermion(L_he, Btier, Ftier, _Hsys.dimensions, Nado, sup_dim, parity, Bbath, Fbath, hierarchy)
end

_getBtier(M::M_Boson_Fermion) = M.Btier
_getFtier(M::M_Boson_Fermion) = M.Ftier

export AbstractHEOMLSMatrix

@doc raw"""
    abstract type AbstractHEOMLSMatrix{T}

Abstract type for all HEOM liouvillian superoperators
"""
abstract type AbstractHEOMLSMatrix{T} end

function Base.getproperty(M::AbstractHEOMLSMatrix, key::Symbol)
    # a comment here to avoid bad render by JuliaFormatter
    if key === :dims
        return dimensions_to_dims(getfield(M, :dimensions))
    else
        return getfield(M, key)
    end
end

@doc raw"""
    (M::AbstractHEOMLSMatrix)(p, t)

Calculate the time-dependent AbstractHEOMLSMatrix at time `t` with parameters `p`.

# Arguments
- `p`: The parameters of the time-dependent coefficients.
- `t`: The time at which the coefficients are evaluated.

# Returns
- The output HEOMLS matrix.
"""
function (M::AbstractHEOMLSMatrix)(p, t)
    # We put 0 in the place of `u` because the time-dependence doesn't depend on the state
    update_coefficients!(M.data, 0, p, t)
    return _reset_HEOMLS_data(M, MatrixOperator(concretize(M.data)))
end

(M::AbstractHEOMLSMatrix)(t) = M(nothing, t)

QuantumToolbox._float_type(M::AbstractHEOMLSMatrix) = _float_type(eltype(M))
QuantumToolbox._complex_float_type(M::AbstractHEOMLSMatrix) = _complex_float_type(eltype(M))

_get_SciML_matrix_wrapper(M::AbstractArray) = QuantumToolbox.get_typename_wrapper(M){eltype(M)}
_get_SciML_matrix_wrapper(M::MatrixOperator) = _get_SciML_matrix_wrapper(M.A)
_get_SciML_matrix_wrapper(M::ScaledOperator) = _get_SciML_matrix_wrapper(M.L)
_get_SciML_matrix_wrapper(M::AddedOperator) = _get_SciML_matrix_wrapper(M.ops[1])
_get_SciML_matrix_wrapper(M::TensorProductOperator) = _get_SciML_matrix_wrapper(M.ops[2].A) # take the system superoperator part for the reference
_get_SciML_matrix_wrapper(M::AbstractHEOMLSMatrix) = _get_SciML_matrix_wrapper(M.data)

# equal to : sparse(vec(system_identity_matrix))
function _Tr(T::Type{<:Number}, dimensions::Dimensions, N::Int)
    D = prod(dimensions)
    return SparseVector(N * D^2, [1 + n * (D + 1) for n in 0:(D-1)], ones(T, D))
end
_Tr(M::AbstractHEOMLSMatrix) = _Tr(eltype(M), M.dimensions, M.N)

function HandleMatrixType(
    M::AbstractQuantumObject,
    MatrixName::String = "";
    type::T = nothing,
) where {T<:Union{Nothing,Operator,SuperOperator}}
    if (type isa Nothing)
        (M.type isa Operator) ||
            (M.type isa SuperOperator) ||
            error("The matrix $(MatrixName) should be either an Operator or a SuperOperator.")
    elseif M.type != type
        error("The matrix $(MatrixName) should be an $(type).")
    end
    return M
end
function HandleMatrixType(
    M::AbstractQuantumObject,
    dimensions::Dimensions,
    MatrixName::String = "";
    type::T = nothing,
) where {T<:Union{Nothing,Operator,SuperOperator}}
    if M.dimensions == dimensions
        return HandleMatrixType(M, MatrixName; type = type)
    else
        error("The dimensions of $(MatrixName) should be: $(_get_dims_string(dimensions))")
    end
end
HandleMatrixType(
    M,
    dimensions::Dimensions,
    MatrixName::String = "";
    type::T = nothing,
) where {T<:Union{Nothing,Operator,SuperOperator}} = HandleMatrixType(M, MatrixName; type = type)
HandleMatrixType(M, MatrixName::String = ""; type::T = nothing) where {T<:Union{Nothing,Operator,SuperOperator}} =
    error("HierarchicalEOM doesn't support matrix $(MatrixName) with type : $(typeof(M))")

# change the type of `ADOs` to match the type of HEOMLS matrix
_HandleVectorType(M::AbstractHEOMLSMatrix, V::SparseVector) = _HandleVectorType(_get_SciML_matrix_wrapper(M), V)
_HandleVectorType(M::Type{<:SparseMatrixCSC}, V::SparseVector) = Vector{_complex_float_type(eltype(M))}(V)

_HandleTraceVectorType(M::AbstractHEOMLSMatrix, V::SparseVector) =
    _HandleTraceVectorType(_get_SciML_matrix_wrapper(M), V)
_HandleTraceVectorType(M::Type{<:SparseMatrixCSC}, V::SparseVector) = V

_HandleSteadyStateMatrix(M::AbstractHEOMLSMatrix{<:MatrixOperator{T,MT}}, b::AbstractVector{T}) where {T<:Number,MT<:SparseMatrixCSC} =
    M.data.A + _SteadyStateConstraint(T, prod(M.dimensions), size(M, 1))
_HandleSteadyStateMatrix(M::AbstractHEOMLSMatrix{<:AbstractSciMLOperator{T}}, b::AbstractVector{T}) where {T<:Number} =
    cache_operator(M.data + _SteadyStateConstraint(eltype(M), prod(M.dimensions), size(M, 1)), b)

# this adds the trace == 1 constraint for reduced density operator during linear solve of steadystate
_SteadyStateConstraint(T::Type{<:Number}, D::Int, S::Int) =
    sparse(ones(T, D), [(n - 1) * (D + 1) + 1 for n in 1:D], ones(T, D), S, S)

function _check_sys_dim_and_ADOs_num(A, B)
    if (A.dimensions != B.dimensions)
        error("Inconsistent system dimensions.")
    end

    if (A.N != B.N)
        error("Inconsistent number of ADOs (\"N\").")
    end
end

_check_parity(A, B) = (typeof(A.parity) != typeof(B.parity)) ? error("Inconsistent parity.") : nothing

function _get_pkg_version(pkg_name::String)
    D = Pkg.dependencies()
    for uuid in keys(D)
        if D[uuid].name == pkg_name
            return D[uuid].version
        end
    end
end

@doc raw"""
    HierarchicalEOM.print_logo(io::IO=stdout)

Print the Logo of HierarchicalEOM package
"""
function print_logo(io::IO = stdout)
    default = :default
    green = 28
    blue = 63
    purple = 133
    red = 124

    printstyled(io, "                                   __", color = green, bold = true)
    printstyled(io, "\n")

    printstyled(io, "                                  /  \\", color = green, bold = true)
    printstyled(io, "\n")

    printstyled(io, " __     __                     ", color = default)
    printstyled(io, "__", color = red, bold = true)
    printstyled(io, " \\__/ ", color = green, bold = true)
    printstyled(io, "__", color = purple, bold = true)
    printstyled(io, "\n")

    printstyled(io, "|  |   |  |                   ", color = default)
    printstyled(io, "/  \\", color = red, bold = true)
    printstyled(io, "    ", color = default)
    printstyled(io, "/  \\", color = purple, bold = true)
    printstyled(io, "\n")

    printstyled(io, "|  |   |  | ______   ______   ", color = default)
    printstyled(io, "\\__/", color = red, bold = true)
    printstyled(io, "_  _", color = default)
    printstyled(io, "\\__/", color = purple, bold = true)
    printstyled(io, "\n")

    printstyled(io, "|  |___|  |/  __  \\ /  ", color = default)
    printstyled(io, "__", color = blue, bold = true)
    printstyled(io, "  \\ / '   \\/     \\", color = default)
    printstyled(io, "\n")

    printstyled(io, "|   ___   |  |__)  |  ", color = default)
    printstyled(io, "/  \\", color = blue, bold = true)
    printstyled(io, "  |    _     _   |", color = default)
    printstyled(io, "\n")

    printstyled(io, "|  |   |  |   ____/| ", color = default)
    printstyled(io, "(    )", color = blue, bold = true)
    printstyled(io, " |   / \\   / \\  |", color = default)
    printstyled(io, "\n")

    printstyled(io, "|  |   |  |  |____ |  ", color = default)
    printstyled(io, "\\__/", color = blue, bold = true)
    printstyled(io, "  |  |   | |   | |", color = default)
    printstyled(io, "\n")

    printstyled(io, "|__|   |__|\\______) \\______/|__|   |_|   |_|", color = default)
    return printstyled(io, "\n")
end

@doc raw"""
    HierarchicalEOM.versioninfo(io::IO=stdout)

Command line output of information on HierarchicalEOM, dependencies, and system information, same as [`HierarchicalEOM.about`](@ref).
"""
function versioninfo(io::IO = stdout)
    cpu = Sys.cpu_info()
    BLAS_info = BLAS.get_config().loaded_libs[1]
    Sys.iswindows() ? OS_name = "Windows" : Sys.isapple() ? OS_name = "macOS" : OS_name = Sys.KERNEL

    # print the logo of HEOM package
    print("\n")
    print_logo(io)

    # print introduction
    println(
        io,
        "\n",
        "Julia framework for Hierarchical Equations of Motion\n",
        "≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡\n",
        "Copyright © QuTiP team 2023 and later.\n",
        "Lead  developer : Yi-Te Huang\n",
        "Other developers:\n",
        "    Simon Cross, Neill Lambert, Po-Chen Kuo and Shen-Liang Yang\n",
    )

    # print package information
    println(
        io,
        "Package information:\n",
        "====================================\n",
        "Julia              Ver. $(VERSION)\n",
        "HierarchicalEOM    Ver. $(_get_pkg_version("HierarchicalEOM"))\n",
        "QuantumToolbox     Ver. $(_get_pkg_version("QuantumToolbox"))\n",
        "SciMLOperators     Ver. $(_get_pkg_version("SciMLOperators"))\n",
        "LinearSolve        Ver. $(_get_pkg_version("LinearSolve"))\n",
        "OrdinaryDiffEqCore Ver. $(_get_pkg_version("OrdinaryDiffEqCore"))\n",
    )

    # print System information
    println(
        io,
        "System information:\n",
        "====================================\n",
        """OS       : $(OS_name) ($(Sys.MACHINE))\n""",
        """CPU      : $(length(cpu)) × $(cpu[1].model)\n""",
        """Memory   : $(round(Sys.total_memory() / 2 ^ 30, digits=3)) GB\n""",
        """WORD_SIZE: $(Sys.WORD_SIZE)\n""",
        """LIBM     : $(Base.libm_name)\n""",
        """LLVM     : libLLVM-$(Base.libllvm_version) ($(Sys.JIT), $(Sys.CPU_NAME))\n""",
        """BLAS     : $(basename(BLAS_info.libname)) ($(BLAS_info.interface))\n""",
        """Threads  : $(Threads.nthreads()) (on $(Sys.CPU_THREADS) virtual cores)\n""",
    )

    # print citation information
    println(
        io,
        "+----------------------------------------------------+\n",
        "| Please cite HierarchicalEOM.jl in your publication |\n",
        "+----------------------------------------------------+\n",
        "For your convenience, a bibtex reference can be easily generated using `HierarchicalEOM.cite()`.\n",
    )
    return nothing
end

@doc raw"""
    HierarchicalEOM.about(io::IO=stdout)

Command line output of information on HierarchicalEOM, dependencies, and system information, same as [`HierarchicalEOM.versioninfo`](@ref).
"""
about(io::IO = stdout) = versioninfo(io)

@doc raw"""
    HierarchicalEOM.cite(io::IO = stdout)

Command line output of citation information and bibtex generator for `HierarchicalEOM.jl`.
"""
function cite(io::IO = stdout)
    citation = raw"""
    @article{HierarchicalEOM.jl2023,
      title = {An efficient {J}ulia framework for hierarchical equations of motion in open quantum systems},
      author = {Huang, Yi-Te and Kuo, Po-Chen and Lambert, Neill and Cirio, Mauro and Cross, Simon and Yang, Shen-Liang and Nori, Franco and Chen, Yueh-Nan},
      journal = {Communications Physics},
      publisher = {Nature Portfolio},
      volume = {6},
      number = {1},
      pages = {313},
      month = {Oct},
      year = {2023},
      doi = {10.1038/s42005-023-01427-2},
      url = {https://doi.org/10.1038/s42005-023-01427-2}
    }
    """
    print(io, citation)
    QuantumToolbox.cite(io)
    return nothing
end

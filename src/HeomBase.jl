abstract type AbstractHEOMLSMatrix{T} end

_FType(::AbstractHEOMLSMatrix{<:AbstractArray{T}}) where {T<:Number} = _FType(T)
_CType(::AbstractHEOMLSMatrix{<:AbstractArray{T}}) where {T<:Number} = _CType(T)

# equal to : sparse(vec(system_identity_matrix))
function _Tr(dims::SVector, N::Int)
    D = prod(dims)
    return SparseVector(N * D^2, [1 + n * (D + 1) for n in 0:(D-1)], ones(ComplexF64, D))
end
function _Tr(M::AbstractHEOMLSMatrix{T}) where {T<:SparseMatrixCSC}
    D = prod(M.dims)
    return SparseVector(M.N * D^2, [1 + n * (D + 1) for n in 0:(D-1)], ones(eltype(M), D))
end

function HandleMatrixType(M::QuantumObject, MatrixName::String = ""; type::QuantumObjectType = Operator)
    if M.type == type
        return M
    else
        error("The matrix $(MatrixName) should be an $(type).")
    end
end
function HandleMatrixType(M::QuantumObject, dims::SVector, MatrixName::String = ""; type::QuantumObjectType = Operator)
    if M.dims == dims
        return HandleMatrixType(M, MatrixName; type = type)
    else
        error("The dims of $(MatrixName) should be: $(dims)")
    end
end
HandleMatrixType(M, dims::SVector, MatrixName::String = ""; type::QuantumObjectType = Operator) =
    HandleMatrixType(M, MatrixName; type = type)
HandleMatrixType(M, MatrixName::String = ""; type::QuantumObjectType = Operator) =
    error("HierarchicalEOM doesn't support matrix $(MatrixName) with type : $(typeof(M))")

# change the type of `ADOs` to match the type of HEOMLS matrix
_HandleVectorType(::AbstractHEOMLSMatrix{<:AbstractSparseMatrix{T}}, V::SparseVector) where {T<:Number} =
    Vector{_CType(T)}(V)

function _HandleSteadyStateMatrix(M::AbstractHEOMLSMatrix, S::Int)
    ElType = eltype(M)
    D = prod(M.dims)
    A = copy(M.data)
    A[1, 1:S] .= 0

    # sparse(row_idx, col_idx, values, row_dims, col_dims)
    A += sparse(ones(ElType, D), [(n - 1) * (D + 1) + 1 for n in 1:D], ones(ElType, D), S, S)
    return A
end

function _check_sys_dim_and_ADOs_num(A, B)
    if (A.dims != B.dims)
        error("Inconsistent system dimension (\"dims\").")
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

"""
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

"""
    HierarchicalEOM.versioninfo(io::IO=stdout)

Command line output of information on HierarchicalEOM, dependencies, and system informations.
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
        "LinearSolve        Ver. $(_get_pkg_version("LinearSolve"))\n",
        "OrdinaryDiffEqCore Ver. $(_get_pkg_version("OrdinaryDiffEqCore"))\n",
    )

    # print System information
    println(io, "System information:")
    println(io, "====================================")
    println(io, """OS       : $(OS_name) ($(Sys.MACHINE))""")
    println(io, """CPU      : $(length(cpu)) × $(cpu[1].model)""")
    println(io, """Memory   : $(round(Sys.total_memory() / 2 ^ 30, digits=3)) GB""")
    println(io, """WORD_SIZE: $(Sys.WORD_SIZE)""")
    println(io, """LIBM     : $(Base.libm_name)""")
    println(io, """LLVM     : libLLVM-$(Base.libllvm_version) ($(Sys.JIT), $(Sys.CPU_NAME))""")
    println(io, """BLAS     : $(basename(BLAS_info.libname)) ($(BLAS_info.interface))""")
    println(io, """Threads  : $(Threads.nthreads()) (on $(Sys.CPU_THREADS) virtual cores)""")
    return print(io, "\n")
end

const PROGBAR_OPTIONS = Dict(:barlen=>20, :color=>:green, :showspeed=>true)

function isValidMatrixType(M, dim::Int=0)
    if typeof(M) <: AbstractMatrix 
        if dim > 0
            if size(M) == (dim, dim)
                return true
            else
                @warn "The size of input matrix should be: ($(M.dim), $(M.dim))."
                return false
            end
        elseif dim == 0
            N1, N2 = size(M)
            if N1 == N2
                return true
            else
                @warn "The size of input matrix should be squared matrix."
                return false
            end
        else
            error("The \"dim\" is not correct.")
        end

    else
        @warn "Heom doesn't support matrix type : $(typeof(M))"
        return false
    end
end

function _get_pkg_version(pkg_name::String)
    D = Pkg.dependencies()
    for uuid in keys(D)
        if D[uuid].name == pkg_name
            return D[uuid].version
        end
    end
end

"""
    Heom.versioninfo(io::IO=stdout)
Command line output of information on Heom, dependencies, and system informations.
"""
function versioninfo(io::IO=stdout)
    cpu = Sys.cpu_info()
    BLAS_config = BLAS.get_config()
    BLAS_build_flags = " (" * join(string.(BLAS_config.build_flags), ", ") * ")"

    # print introduction
    println(io,
        "\n",
        "Heom.jl:\n",
        "Julia framework for Hierarchical Equations of Motion\n",
        "≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡\n",
        "Copyright © NCKU-QFORT 2022 and later.\n",
        "Lead  developer : Yi-Te Huang\n",
        "Other developers: Po-Chen Kuo and Shen-Liang Yang\n"
    )

    # print package informations
    println(io,
        "Package informations:\n",
        "===================================\n",
        "Heom              Version: 2.0.0\n",
        "LinearSolve       Version: $(_get_pkg_version("LinearSolve"))\n",
        "OrdinaryDiffEq    Version: $(_get_pkg_version("OrdinaryDiffEq"))\n",
        "ProgressMeter     Version: $(_get_pkg_version("ProgressMeter"))\n",
        "DistributedArrays Version: $(_get_pkg_version("DistributedArrays"))\n"
    )

    # print System informations
    println(io,
        "System informations:\n",
        "===================================\n",
        "Julia Version: $(VERSION)"
    )
    println(io, 
        "OS       : ", Sys.iswindows() ? "Windows" : Sys.isapple() ? "macOS" : Sys.KERNEL, " (", Sys.MACHINE, ")"
    )
    println(io, 
        "CPU      : ", length(cpu), " × ", cpu[1].model
    )
    println(io,
        "Memory   : ", "$(round(Sys.total_memory() / 2 ^ 30, digits=3)) GB"
    )
    println(io, 
        "WORD_SIZE: ", Sys.WORD_SIZE
    )
    println(io, 
        "LIBM     : ", Base.libm_name
    )
    println(io, 
        "LLVM     : ", "libLLVM-", Base.libllvm_version, " (", Sys.JIT, ", ", Sys.CPU_NAME, ")"
    )
    println(io,
        "BLAS     : ", BLAS.libblastrampoline, BLAS_build_flags
    )
end
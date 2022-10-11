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
    Heom.print_logo(io::IO=stdout)
Print the Logo of Heom package
"""
function print_logo(io::IO=stdout)
    default = Crayon(foreground = :default)
    red     = Crayon(foreground = (203,  60,  51))
    green   = Crayon(foreground = ( 56, 152,  38))
    blue    = Crayon(foreground = ( 64,  99, 216))
    purple  = Crayon(foreground = (149,  88, 178))
    
    print(io, green, "                                   __")
    print(io, "\n")
    
    print(io, green, "                                  /  \\")
    print(io, "\n")
    
    print(io, default, " __     __                     ")
    print(io, red, "__")
    print(io, green, " \\__/ ")
    print(io, purple, "__")
    print(io, "\n")
    
    print(io, default, "|  |   |  |                   ")
    print(io, red, "/  \\")
    print(io, default, "    ")
    print(io, purple, "/  \\")
    print(io, "\n")
    
    print(io, default, "|  |   |  | ______   ______   ")
    print(io, red, "\\__/")
    print(io, default, "_  _")
    print(io, purple, "\\__/")
    print(io, "\n")

    print(io, default, "|  |___|  |/  __  \\ /  ")
    print(io, blue, "__")
    print(io, default, "  \\ / '   \\/     \\")
    print(io, "\n")

    print(io, default, "|   ___   |  |__)  |  ")
    print(io, blue, "/  \\")
    print(io, default, "  |    _     _   |")
    print(io, "\n")

    print(io, default, "|  |   |  |   ____/| ")
    print(io, blue, "(    )")
    print(io, default, " |   / \\   / \\  |")
    print(io, "\n")

    print(io, default, "|  |   |  |  |____ |  ")
    print(io, blue, "\\__/")
    print(io, default, "  |  |   | |   | |")
    print(io, "\n")

    print(io, default, "|__|   |__|\\______) \\______/|__|   |_|   |_|")
    print(io, "\n")
end

"""
    Heom.versioninfo(io::IO=stdout)
Command line output of information on Heom, dependencies, and system informations.
"""
function versioninfo(io::IO=stdout)
    cpu = Sys.cpu_info()
    BLAS_config = BLAS.get_config()
    BLAS_build_flags = " (" * join(string.(BLAS_config.build_flags), ", ") * ")"

    # print the logo of Heom package
    print("\n")
    print_logo(io)

    # print introduction
    println(io,
        "\n",
        "Julia framework for Hierarchical Equations of Motion\n",
        "≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡\n",
        "Copyright © NCKU-QFORT 2022 and later.\n",
        "Lead  developer : Yi-Te Huang\n",
        "Other developers:\n",
        "    Neill Lambert, Po-Chen Kuo and Shen-Liang Yang\n"
    )

    # print package informations
    println(io,
        "Package informations:\n",
        "===================================\n",
        "Heom              Version: 0.1.0\n",
        "DistributedArrays Version: $(_get_pkg_version("DistributedArrays"))\n",
        "LinearSolve       Version: $(_get_pkg_version("LinearSolve"))\n",
        "OrdinaryDiffEq    Version: $(_get_pkg_version("OrdinaryDiffEq"))\n",
        "ProgressMeter     Version: $(_get_pkg_version("ProgressMeter"))\n",
        "SnoopPrecompile   Version: $(_get_pkg_version("SnoopPrecompile"))"
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
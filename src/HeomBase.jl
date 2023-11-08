const PROGBAR_OPTIONS = Dict(:barlen=>20, :color=>:green, :showspeed=>true)

function HandleMatrixType(M, dim::Int=0, MatrixName::String="")
    error("HierarchicalEOM doesn't support matrix $(MatrixName) with type : $(typeof(M))")
end

function HandleMatrixType(M::AbstractMatrix, dim::Int=0, MatrixName::String="")
    if dim > 0
        if size(M) == (dim, dim)
            return copy(M)
        else
            error("The size of matrix $(MatrixName) should be: ($(dim), $(dim)).")
        end
    elseif dim == 0
        N1, N2 = size(M)
        if N1 == N2
            return copy(M)
        else
            error("The size of matrix $(MatrixName) should be squared matrix.")
        end
    end
end

function _HandleFloatType(ElType::Type{T}, V::StepRangeLen) where T <: Number
    if real(ElType) == Float32
        return StepRangeLen(Float32(V.ref), Float32(V.step), Int32(V.len), Int64(V.offset))
    else
        return StepRangeLen(Float64(V.ref), Float64(V.step), Int64(V.len), Int64(V.offset))
    end
end

function _HandleFloatType(ElType::Type{T}, V::Any) where T <: Number
    FType = real(ElType)
    if eltype(V) == FType
        return V
    else
        convert.(FType, V)
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
    HierarchicalEOM.print_logo(io::IO=stdout)
Print the Logo of HierarchicalEOM package
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
    HierarchicalEOM.versioninfo(io::IO=stdout)
Command line output of information on HierarchicalEOM, dependencies, and system informations.
"""
function versioninfo(io::IO=stdout)
    cpu = Sys.cpu_info()
    BLAS_info = BLAS.get_config().loaded_libs[1]

    # print the logo of HEOM package
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
        "    Simon Cross, Neill Lambert, Po-Chen Kuo and Shen-Liang Yang\n"
    )

    # print package informations
    println(io,
        "Package information:\n",
        "====================================\n",
        "HierarchicalEOM   Ver. $(_get_pkg_version("HierarchicalEOM"))\n",
        "LinearSolve       Ver. $(_get_pkg_version("LinearSolve"))\n",
        "OrdinaryDiffEq    Ver. $(_get_pkg_version("OrdinaryDiffEq"))\n",
        "SteadyStateDiffEq Ver. $(_get_pkg_version("SteadyStateDiffEq"))\n",
        "SciMLOperators    Ver. $(_get_pkg_version("SciMLOperators"))\n",
        "FastExpm          Ver. $(_get_pkg_version("FastExpm"))\n",
        "JLD2              Ver. $(_get_pkg_version("JLD2"))\n",
        "ProgressMeter     Ver. $(_get_pkg_version("ProgressMeter"))\n",
        "PrecompileTools   Ver. $(_get_pkg_version("PrecompileTools"))\n"        
    )

    # print System informations
    println(io,
        "System information:\n",
        "====================================\n",
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
        "BLAS     : ", basename(BLAS_info.libname), " (", BLAS_info.interface, ")"
    )
    print(io, "\n")
end
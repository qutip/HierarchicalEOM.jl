# julia --project=@. make.jl
push!(LOAD_PATH, "../src/")
using Documenter
import Literate
using HierarchicalEOM

const DRAFT = false # set `true` to disable cell evaluation

doc_output_path = abspath(joinpath(@__DIR__, "src", "examples"))
if !isdir(doc_output_path) mkdir(doc_output_path) end

# Generate page: Quick Start
QS_source_file = abspath(joinpath(@__DIR__, "..",  "examples", "quick_start.jl"))
Literate.markdown(QS_source_file, doc_output_path)

const PAGES = Any[
    "Home" => Any[
        "Introduction" => "index.md",
        "Installation" => "install.md",
        "Quick Start"  => "examples/quick_start.md",
        "Cite HierarchicalEOM" => "cite.md"
    ],
    "Bosonic Bath" => "bosonic_bath.md",
    "Fermionic Bath" => "fermionic_bath.md",
    "Auxiliary Density Operators" => "ADOs.md",
    "HEOMLS Matrices" => Any[
        "Introduction" => "heom_matrix/intro.md",
        "HEOMLS for Schrödinger Equation" => "heom_matrix/schrodinger_eq.md",
        "HEOMLS for Bosonic Bath" => "heom_matrix/M_Boson.md",
        "HEOMLS for Fermionic Bath" => "heom_matrix/M_Fermion.md",
        "HEOMLS for Bosonic and Fermionic Bath" => "heom_matrix/M_Boson_Fermion.md",
        "HEOMLS for Master Equation" => "heom_matrix/master_eq.md",
    ],
    "Hierarchy Dictionary" => "hierarchy_dictionary.md",
    "Time Evolution" => "time_evolution.md",
    "Stationary State" => "stationary_state.md",
    "Spectrum" => "spectrum.md",
    # "Examples" => Any[],
    "Library" => Any[
        "HEOM API" => "lib/heom_api.md",
        "Bath" => "lib/bath.md",
        "Bath Correlation Functions" => "lib/corr_func.md",
        "Physical Analysis Functions" => "lib/phys_analysis.md",
        "Misc." => "lib/misc.md"
    ]
]

makedocs(
    sitename = "Documentation | HierarchicalEOM.jl",
    pages  = PAGES,
    format = Documenter.HTML(
        prettyurls = (get(ENV, "CI", nothing) == "true"),    
        ansicolor  = true
    ),
    draft  = DRAFT
)

deploydocs(
    repo="github.com/NCKU-QFort/HierarchicalEOM.jl.git",
    devbranch = "main"
)